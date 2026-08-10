#ifndef Pipe1_PIDM4_cxx
#define Pipe1_PIDM4_cxx

#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActParticle.h"
#include "ActSRIM.h"
#include "ActSilData.h"
#include "ActSilMatrix.h"
#include "ActSilSpecs.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TH2.h"
#include "TH2D.h"
#include "TMath.h"
#include "TString.h"

#include <search.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../../../PrettyStyle.C"
#include "../Utils.h"

void Pipe1_PIDM4(const std::string& beam, const std::string& target, const std::string& light)
{
    // Get file from pipe0
    TString infile = TString::Format("./Outputs/SelectorM4_%s.root", beam.c_str());
    // ROOT::EnableImplicitMT();
    ROOT::RDataFrame df("Final_Tree", infile.Data());

    // Get silicon matrix and silspecs to get particle to go to lateral silicons
    std::string filenameSMleft {"../SilVetos/Outputs/Dists/sms_l0.root"};
    auto fileSMleft {new TFile {filenameSMleft.c_str()}};
    ActPhysics::SilMatrix* smleft =
        fileSMleft->Get<ActPhysics::SilMatrix>("sm5"); // matrix for good distance of left wall
    double silCentreLeft = smleft->GetMeanZ({4, 5});
    double beamOffsetLeft {3.36}; // mm offset of beam with respect to sils 4 and 5 off left wall (need to lower beam
                                  // that amount respect of silicons)
    const double zVertexMeanLeft {silCentreLeft - beamOffsetLeft}; // beam on left
    std::string filenameSMright {"../SilVetos/Outputs/Dists/sms_r0.root"};
    auto fileSMright {new TFile {filenameSMright.c_str()}};
    ActPhysics::SilMatrix* smright =
        fileSMright->Get<ActPhysics::SilMatrix>("sm5"); // matrix for good distance of right wall
    double silCentreRight = smright->GetMeanZ({4, 5});
    double beamOffsetRight {-2.14};
    const double zVertexMeanRight {silCentreRight - beamOffsetRight}; // beam on right
    std::string filenameSMfront {"../SilVetos/Outputs/Dists/sms_f0.root"};
    auto fileSMfront {new TFile {filenameSMfront.c_str()}};
    ActPhysics::SilMatrix* smfront =
        fileSMfront->Get<ActPhysics::SilMatrix>("sm2"); // matrix for good distance of front wall
    double silCentreFront = smfront->GetMeanZ({7, 6, 4});
    double beamOffsetFront {-2.5};
    const double zVertexMeanFront {silCentreFront - beamOffsetFront}; // beam on front

    auto sils {std::make_shared<ActPhysics::SilSpecs>()};
    sils->ReadFile("../../configs/silspecs.conf");

    ActRoot::InputParser parser {};
    parser.ReadFile("../../configs/detector.conf");
    auto mergerBlock = parser.GetBlock("Merger");
    auto validationZone = mergerBlock->GetDouble("L1ExclusionZone"); // in pad units

    auto driftFactor = mergerBlock->GetDouble("DriftFactor"); // in pad units

    // Read per-run zDrift cuts: columns are run, meanStart, sigmaStart, meanEnd, sigmaEnd
    struct ZDriftCut
    {
        double zMin; // meanStart + nSigma * sigmaStart
        double zMax; // meanEnd   - nSigma * sigmaEnd
    };
    const double nSigmaZDrift {3.0};
    std::map<int, ZDriftCut> zDriftCuts;
    {
        std::ifstream finZ("../L1PID/Outputs/zDrift_cut_perRun.dat");
        if(!finZ.is_open())
            throw std::runtime_error("Could not open zDrift_cut_perRun.dat");
        std::string line;
        while(std::getline(finZ, line))
        {
            if(line.empty() || line[0] == '#')
                continue;
            std::istringstream iss(line);
            int run;
            double meanS, sigmaS, meanE, sigmaE;
            if(iss >> run >> meanS >> sigmaS >> meanE >> sigmaE)
                zDriftCuts[run] = {meanS + nSigmaZDrift * sigmaS, meanE - nSigmaZDrift * sigmaE};
        }
        finZ.close();
    }


    // Now get the index of the particle that hit the silicon to do PID
    auto dfLayer =
        df.Define("StoppedInACTAR",
                  [&](ActRoot::MergerData& m, ActRoot::TPCData& tpc, int LightIdx, float zDrift, float Qave)
                  {
                      // First check to be L1 is needed to stop within the limits of the TPC
                      if(LightIdx < 0)
                          return false;
                      auto cluster = tpc.fClusters[LightIdx]; // copia local, no modifica el original
                      cluster.SortAlongDir(cluster.GetLine().GetDirection());
                      auto& voxels = cluster.GetRefToVoxels();
                      if(voxels.empty())
                          return false;
                      auto stopPoint = voxels.back().GetPosition();
                      bool insideX = (stopPoint.X() > validationZone && stopPoint.X() < (128 - validationZone));
                      bool insideY = (stopPoint.Y() > validationZone && stopPoint.Y() < (128 - validationZone));

                      // Now check if the zDrift and Qave enters the cuts per run
                      bool isInsideZDriftCut = false;
                      bool isInsideQaveCut = false;
                      auto it = zDriftCuts.find(m.fRun);
                      if(it == zDriftCuts.end())
                          isInsideZDriftCut = false; // run not in file → reject
                      else
                          isInsideZDriftCut = zDrift > it->second.zMin && zDrift < it->second.zMax;

                      isInsideQaveCut = Qave > 400;

                      return insideX && insideY && isInsideZDriftCut && isInsideQaveCut;
                  },
                  {"MergerData", "TPCData", "LightIdx", "zDrift", "Qave"})
            .Define(
                "LayerAndSP",
                [&](ActRoot::TPCData& tpc, ActRoot::SilData& silIn,
                    int LightIdx) -> std::pair<std::string, ROOT::Math::XYZPointF>
                {
                    std::pair<std::string, ROOT::Math::XYZPointF> res {"", {}};
                    if(LightIdx < 0)
                        return res;

                    ActRoot::SilData sil {silIn};
                    sil.ApplyFinerThresholds(sils);

                    auto& cluster = tpc.fClusters[LightIdx];
                    auto line = cluster.GetLine();
                    line.Scale(Utils::scaleXY, Utils::scaleZ);

                    double minDistance {1e8};
                    const double kSilMatchTol {
                        70.}; // mm -- distance between front and sides is 71 mm. So has to be lower than that

                    for(const std::string layer : {"f0", "l0", "r0"})
                    {
                        auto it = sil.fSiN.find(layer);
                        if(it == sil.fSiN.end() || it->second.empty())
                            continue; // no hit en esta layer

                        auto pointLayer {sils->GetLayer(layer).GetPoint()};
                        auto silIdx {it->second[0]};
                        auto silHitPosition {sils->GetLayer(layer).GetPlacements().at(silIdx)};

                        ROOT::Math::XYZPointF pseudoSP {};
                        if(layer == "l0")
                            pseudoSP = {static_cast<float>(silHitPosition.first), static_cast<float>(pointLayer.Y()),
                                        static_cast<float>(smleft->GetMeanZ({silIdx}) - beamOffsetLeft)};
                        else if(layer == "r0")
                            pseudoSP = {static_cast<float>(silHitPosition.first), static_cast<float>(pointLayer.Y()),
                                        static_cast<float>(smright->GetMeanZ({silIdx}) - beamOffsetRight)};
                        else if(layer == "f0")
                            pseudoSP = {static_cast<float>(silHitPosition.first), static_cast<float>(pointLayer.Y()),
                                        static_cast<float>(smfront->GetMeanZ({silIdx}) - beamOffsetFront)};
                        // if(layer == "l0")
                        //     pseudoSP = {static_cast<float>(pointLayer.X()), static_cast<float>(pointLayer.Y()), 0};
                        // else if(layer == "r0")
                        //     pseudoSP = {static_cast<float>(pointLayer.X()), static_cast<float>(pointLayer.Y()), 0};
                        // else if(layer == "f0")
                        //     pseudoSP = {static_cast<float>(pointLayer.X()), static_cast<float>(pointLayer.Y()), 0};

                        ROOT::Math::XYZPoint posibleSP = {0, 0, 0};
                        float distance {};
                        // if(layer == "l0" || layer == "r0")
                        //     distance = (posibleSP - pseudoSP).R();
                        if(layer == "l0" || layer == "r0")
                        {
                            posibleSP = line.MoveToY(pointLayer.Y());
                            // posibleSP.SetZ(0); // We have a bad matrix for the front silicon, so just check in XY
                            // plane
                            distance = (posibleSP - pseudoSP).R();
                        }
                        else if(layer == "f0")
                        {
                            posibleSP = line.MoveToX(pointLayer.X());
                            // posibleSP.SetZ(0); // We have a bad matrix for the front silicon, so just check in XY
                            // plane
                            distance = (posibleSP - pseudoSP).R();
                        }

                        if(distance < kSilMatchTol && distance < minDistance)
                        {
                            minDistance = distance;
                            res = {layer, pseudoSP};
                        }
                    }
                    return res;
                },
                {"TPCData", "SilData", "LightIdx"})
            .Define("SilLayer", [](const std::pair<std::string, ROOT::Math::XYZPointF>& p) { return p.first; },
                    {"LayerAndSP"})
            .Define("pseudoSP", [](const std::pair<std::string, ROOT::Math::XYZPointF>& p) { return p.second; },
                    {"LayerAndSP"})
            .Define("Layer",
                    [](const std::string& silLayer, bool stopped)
                    {
                        // Si la particula cumple las condiciones de L1 (parada validada dentro de ACTAR),
                        // se queda como L1 aunque tambien exista un matching geometrico con algun silicio.
                        if(stopped)
                            return std::string("L1");
                        if(silLayer != "")
                            return silLayer;    // si no es L1, el matching de silicio manda, si existe
                        return std::string(""); // ni silicio ni parada valida -> se descarta
                    },
                    {"SilLayer", "StoppedInACTAR"});

    auto dfWithIndex =
        dfLayer
            .Filter([](const std::string& layer) { return layer != ""; },
                    {"Layer"}) // descarta si no matchea ninguna layer
            .Define("TrackLength",
                    [&](ActRoot::TPCData& tpc, int LightIdx, ROOT::Math::XYZPointF pseudoSP, std::string layer)
                    {
                        auto& cluster = tpc.fClusters[LightIdx];
                        auto line = cluster.GetLine();
                        line.Scale(Utils::scaleXY, Utils::scaleZ);
                        auto rp = tpc.fRPs.front();
                        Utils::ScalePoint(rp); // Scale RP to same units as line and pseudoSP
                        if(layer == "L1")
                        {
                            auto& voxels = cluster.GetRefToVoxels();
                            auto stopPoint = voxels.back().GetPosition();
                            Utils::ScalePoint(stopPoint); // misma escala que se usa en Qave, para consistencia
                            return static_cast<float>((stopPoint - rp).R());
                        }
                        else
                        {
                            auto sp = line.MoveToY(pseudoSP.Y());
                            return static_cast<float>((sp - rp).R());
                        }
                    },
                    {"TPCData", "LightIdx", "pseudoSP", "Layer"})
            .Define("TrackLengthPadUnits",
                    [&](ActRoot::TPCData& tpc, int LightIdx, ROOT::Math::XYZPointF pseudoSP, std::string layer)
                    {
                        auto& cluster = tpc.fClusters[LightIdx];
                        auto line = cluster.GetLine();
                        auto rp = tpc.fRPs.front();
                        if(layer == "L1")
                        {
                            auto& voxels = cluster.GetRefToVoxels();
                            auto stopPoint = voxels.back().GetPosition();
                            return static_cast<float>((stopPoint - rp).R());
                        }
                        else
                        {
                            auto sp = line.MoveToY(pseudoSP.Y());
                            sp.SetXYZ(sp.X() / 2, sp.Y() / 2, sp.Z() / driftFactor); // back to pad units
                            return static_cast<float>((sp - rp).R());
                        }
                    },
                    {"TPCData", "LightIdx", "pseudoSP", "Layer"})
            .Define("BeamIdx",
                    [](ActRoot::TPCData& tpc)
                    {
                        int beamIdx = -1;
                        for(int i = 0; i < tpc.fClusters.size(); ++i)
                        {
                            if(tpc.fClusters[i].GetIsBeamLike())
                            {
                                beamIdx = i;
                                break;
                            }
                        }
                        return beamIdx;
                    },
                    {"TPCData"});

    auto hIdx =
        dfWithIndex.Histo1D({"hIdx", "Idx of particle hitting lateral silicon;Index;Counts", 7, -2, 5}, "LightIdx");

    // Now define Qave for that index to do PID plot
    auto dfWithQave =
        dfWithIndex
            .Define("SilE",
                    [&](ActRoot::SilData& sil, ActRoot::ModularData& mod, std::string layer)
                    {
                        if(layer == "L1")
                            return -1.f; // no hay energia de silicio para eventos parados dentro
                        sil.ApplyFinerThresholds(sils);
                        return sil.fSiE.at(layer)[0];
                    },
                    {"SilData", "ModularData", "Layer"})
            .Define("ThetaLab",
                    [&](ActRoot::TPCData& tpc, int Lightidx, int BeamIdx)
                    {
                        auto cluster = tpc.fClusters[Lightidx];
                        auto beam = tpc.fClusters[BeamIdx];
                        auto theta = Utils::GetTheta3D(beam.GetLine().GetDirection(), cluster.GetLine().GetDirection());
                        return theta; // in degrees
                    },
                    {"TPCData", "LightIdx", "BeamIdx"})
            .Define("Qtot",
                    [&](ActRoot::TPCData& tpc, int LightIdx)
                    {
                        auto& cluster = tpc.fClusters[LightIdx];
                        float Qtot = 0;
                        for(const auto& voxel : cluster.GetRefToVoxels())
                        {
                            Qtot += voxel.GetCharge();
                        }
                        return Qtot;
                    },
                    {"TPCData", "LightIdx"});

    // Plot PID to check everything is ok, first divide in sil layers
    auto dfL0 = dfWithQave.Filter([](const std::string& layer) { return layer == "l0"; }, {"Layer"});
    auto dfR0 = dfWithQave.Filter([](const std::string& layer) { return layer == "r0"; }, {"Layer"});
    auto dfF0 = dfWithQave.Filter([](const std::string& layer) { return layer == "f0"; }, {"Layer"});
    auto dfL1 = dfWithQave.Filter([](const std::string& layer) { return layer == "L1"; }, {"Layer"});

    auto hPIDl0 =
        dfL0.Histo2D({"hPIDl0", "PID plot Qave vs SilE;Silicon Energy (MeV);Q_{ave} (a.u.)", 100, 0, 30, 100, 0, 500},
                     "SilE", "Qave");
    auto hPIDr0 =
        dfR0.Histo2D({"hPIDr0", "PID plot Qave vs SilE;Silicon Energy (MeV);Q_{ave} (a.u.)", 100, 0, 30, 100, 0, 500},
                     "SilE", "Qave");
    auto hPIDf0 =
        dfF0.Histo2D({"hPIDf0", "PID plot Qave vs SilE;Silicon Energy (MeV);Q_{ave} (a.u.)", 100, 0, 30, 100, 0, 500},
                     "SilE", "Qave");
    auto hPIDl1 = dfL1.Histo2D(
        {"hPIDl1", "PID plot Qtot vs TrackLength;TrackLength (pad);Q_{tot} (a.u.)", 200, 0, 120, 2000, 0, 3e5},
        "TrackLengthPadUnits", "Qtot");

    auto cuts = ActRoot::CutsManager<std::string> {};
    cuts.ReadCut("l0",
                 TString::Format("../../PostAnalysis/Cuts/pid_%s_l0_%s.root", light.c_str(), beam.c_str()).Data());
    cuts.ReadCut("r0",
                 TString::Format("../../PostAnalysis/Cuts/pid_%s_r0_%s.root", light.c_str(), beam.c_str()).Data());
    cuts.ReadCut("f0",
                 TString::Format("../../PostAnalysis/Cuts/pid_%s_f0_%s.root", light.c_str(), beam.c_str()).Data());
    cuts.ReadCut("L1", TString::Format("./Cuts/pid_%s_l1_%s_strict.root", light.c_str(), beam.c_str()).Data());

    std::cout << "Number of events with M=4 and later silicon hit: " << dfWithQave.Count().GetValue() << std::endl;

    // Finally filter the events with the cuts
    auto dfPIDfiltered = dfWithQave.Filter(
        [&cuts](std::string layer, float silE, float Qave, float TL, float Qtot)
        {
            if(layer == "l0")
                return cuts.IsInside("l0", silE, Qave);
            else if(layer == "r0")
                return cuts.IsInside("r0", silE, Qave);
            else if(layer == "L1")
                return cuts.IsInside("L1", TL, Qtot);
            else
                return false;
        },
        {"Layer", "SilE", "Qave", "TrackLengthPadUnits", "Qtot"});

    TCanvas* c0 = new TCanvas("Pipe1_0", "Index of particle hitting lateral silicon", 800, 600);
    hIdx->DrawClone();

    TCanvas* c1 = new TCanvas("Pipe1_1", "PID plot", 1600, 1200);
    c1->Divide(2, 2);
    c1->cd(1);
    hPIDl0->DrawClone("colz");
    cuts.DrawCut("l0");
    c1->cd(2);
    hPIDr0->DrawClone("colz");
    cuts.DrawCut("r0");
    c1->cd(3);
    hPIDf0->DrawClone("colz");
    cuts.DrawCut("f0");
    c1->cd(4);
    hPIDl1->DrawClone("colz");
    cuts.DrawCut("L1");

    // SilLayerMask: bitmask indicando que capas se activaron para el evento
    // bit0=f0, bit1=l0, bit2=r0, bit3=L1 (parada validada dentro de ACTAR)
    auto dfOverlap =
        dfLayer.Define("SilLayerMask",
                       [&](ActRoot::SilData& silIn, bool stoppedInACTAR)
                       {
                           ActRoot::SilData sil {silIn};
                           sil.ApplyFinerThresholds(sils);

                           bool hasF0 = sil.fSiN.count("f0") && !sil.fSiN.at("f0").empty();
                           bool hasL0 = sil.fSiN.count("l0") && !sil.fSiN.at("l0").empty();
                           bool hasR0 = sil.fSiN.count("r0") && !sil.fSiN.at("r0").empty();

                           // bitmask: bit0=f0, bit1=l0, bit2=r0, bit3=L1
                           int mask = (hasF0 ? 1 : 0) | (hasL0 ? 2 : 0) | (hasR0 ? 4 : 0) | (stoppedInACTAR ? 8 : 0);
                           return mask;
                       },
                       {"SilData", "StoppedInACTAR"});
    // std::ofstream out1(TString::Format("./Outputs/Pipe1_PIDM4_F0andSide_butNotSideLayerCount.dat").Data());
    // dfOverlap.Foreach(
    //     [&](ActRoot::MergerData& m, const std::string& layer, int mask)
    //     {
    //         if(layer == "f0" && (mask == 3 || mask == 5)) // f0 and l0 or r0 but not both
    //             m.Stream(out1);
    //     },
    //     {"MergerData", "Layer", "SilLayerMask"});
    // out1.close();
    std::ofstream out1(TString::Format("./Outputs/Pipe1_PIDM4_onlyF0.dat").Data());
    dfOverlap.Foreach(
        [&](ActRoot::MergerData& m, const std::string& layer, int mask)
        {
            if(mask == 1) // f0 and l0 or r0 but not both
                m.Stream(out1);
        },
        {"MergerData", "Layer", "SilLayerMask"});
    out1.close();
    std::ofstream out2(TString::Format("./Outputs/Pipe1_PIDM4_otherEvents.dat").Data());
    dfOverlap.Foreach(
        [&](ActRoot::MergerData& m, const std::string& layer, int mask)
        {
            if(mask == 0) // f0 and l0 or r0 but not both
                m.Stream(out2);
        },
        {"MergerData", "Layer", "SilLayerMask"});
    out2.close();
    std::ofstream out3(TString::Format("./Outputs/Pipe1_PIDM4_L1Events.dat").Data());
    dfOverlap.Foreach(
        [&](ActRoot::MergerData& m, const std::string& layer, int mask)
        {
            if(mask >= 8) // L1
                m.Stream(out3);
        },
        {"MergerData", "Layer", "SilLayerMask"});
    out3.close();

    // Ahora el histograma cubre las 16 combinaciones posibles (f0, l0, r0, L1)
    auto hOverlap =
        dfOverlap.Histo1D({"hOverlap", "Silicon layer co-hits;Category;Counts", 16, -0.5, 15.5}, "SilLayerMask");

    TCanvas* cOverlap = new TCanvas("cOverlap", "f0 vs l0/r0 vs L1 overlap", 1000, 700);
    hOverlap->DrawClone("hist text0");

    // Etiquetas legibles para cada bin (bit0=F0, bit1=L0, bit2=R0, bit3=L1)
    const char* labels[16] = {"None", "F0",    "L0",    "F0+L0",    "R0",    "F0+R0",    "L0+R0",    "F0+L0+R0",
                              "L1",   "F0+L1", "L0+L1", "F0+L0+L1", "R0+L1", "F0+R0+L1", "L0+R0+L1", "F0+L0+R0+L1"};
    for(int i = 0; i < 16; ++i)
        hOverlap->GetXaxis()->SetBinLabel(i + 1, labels[i]);
    hOverlap->DrawClone("hist text0");

    std::cout << "Counts per category:\n";
    for(int i = 1; i <= 16; ++i)
        std::cout << "  " << hOverlap->GetXaxis()->GetBinLabel(i) << ": " << hOverlap->GetBinContent(i) << '\n';

    std::cout << "Raw con hit lateral (l0/r0, con o sin f0): 643 (ya calculado)\n";
    std::cout
        << "Tras geometric matching (Layer != \"\"): "
        << dfWithIndex.Filter([](const std::string& l) { return l == "l0" || l == "r0"; }, {"Layer"}).Count().GetValue()
        << std::endl;
    std::cout << "Tras corte PID: " << dfPIDfiltered.Count().GetValue() << std::endl;

    // std::ofstream goodLateral(TString::Format("./Outputs/Pipe1_PIDM4_goodlateral.dat").Data());
    // dfL0.Foreach([&](ActRoot::MergerData& m) { m.Stream(goodLateral); }, {"MergerData"});
    // dfR0.Foreach([&](ActRoot::MergerData& m) { m.Stream(goodLateral); }, {"MergerData"});

    std::vector<std::pair<int, int>> goodKeys;

    // dfL0.Foreach([&](ActRoot::MergerData& m) { goodKeys.push_back({m.fRun, m.fEntry}); }, {"MergerData"});
    // dfR0.Foreach([&](ActRoot::MergerData& m) { goodKeys.push_back({m.fRun, m.fEntry}); }, {"MergerData"});
    //
    // std::ofstream f0AndLateral_NotInPID(TString::Format("./Outputs/Pipe1_PIDM4_f0AndLateral_NotInPID.dat").Data());
    // dfOverlap.Foreach(
    //     [&](ActRoot::MergerData& m, int mask)
    //     {
    //         if(mask == 2 || mask == 3 || mask == 4 || mask == 5) // l0 or r0 or f0+l0 or f0+r0
    //         {
    //
    //             bool found = false;
    //             int i = 0;
    //             while(i < (int)goodKeys.size())
    //             {
    //                 if(goodKeys[i].first == m.fRun && goodKeys[i].second == m.fEntry)
    //                 {
    //                     found = true;
    //                     break;
    //                 }
    //                 i++;
    //             }
    //
    //             if(!found)
    //             {
    //                 m.Stream(f0AndLateral_NotInPID);
    //             }
    //         }
    //     },
    //     {"MergerData", "SilLayerMask"});
    // f0AndLateral_NotInPID.close();
    //
    // std::ofstream out(TString::Format("./Outputs/Pipe1_PIDM4_side_%s.dat", light.c_str()).Data());
    // dfPIDfiltered.Foreach([&](ActRoot::MergerData& m, const std::string& layer) { m.Stream(out); },
    //                       {"MergerData", "Layer"});
    // out.close();

    // Save dataframe in a .root file
    TString outfile = TString::Format("./Outputs/PIDM4_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str());
    dfPIDfiltered.Snapshot("Final_Tree", outfile);
    std::cout << "Saving Final_Tree in " << outfile << " from Pipe1_PIDM4" << '\n';
}
#endif