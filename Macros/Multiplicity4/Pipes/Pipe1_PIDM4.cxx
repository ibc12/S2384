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

    auto sils {std::make_shared<ActPhysics::SilSpecs>()};
    sils->ReadFile("../../configs/silspecs.conf");


    // Now get the index of the particle that hit the silicon to do PID
    auto dfLayer =
        df.Define("LayerAndSP",
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
                      const double kSilMatchTol {60.}; // mm -- ajustar según resolución angular/geométrica

                      for(const std::string layer : {"l0", "r0", "f0"})
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
                              pseudoSP = {static_cast<float>(pointLayer.X()), static_cast<float>(pointLayer.Y()), 0};

                          auto posibleSP = line.MoveToY(pseudoSP.Y());
                          float distance {};
                          if(layer == "l0" || layer == "r0")
                              distance = (posibleSP - pseudoSP).R();
                          else if(layer == "f0")
                          {
                              posibleSP.SetZ(
                                  0); // We have a bad matrix for the front silicon, so just check in XY plane
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
            .Define("Layer", [](const std::pair<std::string, ROOT::Math::XYZPointF>& p) { return p.first; },
                    {"LayerAndSP"})
            .Define("pseudoSP", [](const std::pair<std::string, ROOT::Math::XYZPointF>& p) { return p.second; },
                    {"LayerAndSP"});
    auto dfWithIndex = dfLayer
                           .Filter([](const std::string& layer) { return layer != ""; },
                                   {"Layer"}) // descarta si no matchea ninguna layer
                           .Define("TrackLength",
                                   [&](ActRoot::TPCData& tpc, int LightIdx, ROOT::Math::XYZPointF pseudoSP)
                                   {
                                       auto& cluster = tpc.fClusters[LightIdx];
                                       auto line = cluster.GetLine();
                                       line.Scale(Utils::scaleXY, Utils::scaleZ);
                                       auto rp = tpc.fRPs.front();
                                       auto sp = line.MoveToY(pseudoSP.Y());
                                       return static_cast<float>((sp - rp).R());
                                   },
                                   {"TPCData", "LightIdx", "pseudoSP"})
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
            .Define("QaveSelectedParticle",
                    [](ActRoot::TPCData& tpc, int idx)
                    {
                        auto cluster = tpc.fClusters[idx];
                        cluster.SortAlongDir(cluster.GetLine().GetDirection());
                        auto voxels = cluster.GetRefToVoxels();
                        float Q = 0.f;
                        for(const auto& v : voxels)
                            Q += v.GetCharge();
                        ROOT::Math::XYZPointF firstPos = voxels.front().GetPosition();
                        ROOT::Math::XYZPointF lastPos = voxels.back().GetPosition();
                        Utils::ScalePoint(firstPos);
                        Utils::ScalePoint(lastPos);
                        float length = (lastPos - firstPos).R();
                        if(length <= 0.0)
                            return -1.f;
                        return Q / length;
                    },
                    {"TPCData", "LightIdx"})
            .Define("SilESelectedParticle",
                    [&](ActRoot::SilData& sil, ActRoot::ModularData& mod, std::string layer)
                    {
                        // std::cout << "------------------------" << std::endl;
                        // std::cout << "Size of hits in layer  before threshold" << layer << ": "
                        //           << sil.fSiE.at(layer).size() << std::endl;
                        sil.ApplyFinerThresholds(sils);
                        // std::cout << "After threshold: " << sil.fSiE.at(layer).size() << std::endl;
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
                    {"TPCData", "LightIdx", "BeamIdx"});

    // Plot PID to check everything is ok, first divide in sil layers
    auto dfL0 = dfWithQave.Filter([](const std::string& layer) { return layer == "l0"; }, {"Layer"});
    auto dfR0 = dfWithQave.Filter([](const std::string& layer) { return layer == "r0"; }, {"Layer"});
    auto dfF0 = dfWithQave.Filter([](const std::string& layer) { return layer == "f0"; }, {"Layer"});

    auto hPIDl0 =
        dfL0.Histo2D({"hPIDl0", "PID plot Qave vs SilE;Silicon Energy (MeV);Q_{ave} (a.u.)", 100, 0, 30, 100, 0, 500},
                     "SilESelectedParticle", "QaveSelectedParticle");
    auto hPIDr0 =
        dfR0.Histo2D({"hPIDr0", "PID plot Qave vs SilE;Silicon Energy (MeV);Q_{ave} (a.u.)", 100, 0, 30, 100, 0, 500},
                     "SilESelectedParticle", "QaveSelectedParticle");
    auto hPIDf0 =
        dfF0.Histo2D({"hPIDf0", "PID plot Qave vs SilE;Silicon Energy (MeV);Q_{ave} (a.u.)", 100, 0, 30, 100, 0, 500},
                     "SilESelectedParticle", "QaveSelectedParticle");

    auto cuts = ActRoot::CutsManager<std::string> {};
    cuts.ReadCut("l0",
                 TString::Format("../../PostAnalysis/Cuts/pid_%s_l0_%s.root", light.c_str(), beam.c_str()).Data());
    cuts.ReadCut("r0",
                 TString::Format("../../PostAnalysis/Cuts/pid_%s_r0_%s.root", light.c_str(), beam.c_str()).Data());
    cuts.ReadCut("f0",
                 TString::Format("../../PostAnalysis/Cuts/pid_%s_f0_%s.root", light.c_str(), beam.c_str()).Data());

    std::cout << "Number of events with M=4 and later silicon hit: " << dfWithQave.Count().GetValue() << std::endl;

    // Finally filter the events with the cuts
    auto dfPIDfiltered = dfWithQave.Filter(
        [&cuts](std::string layer, float silE, float Qave)
        {
            if(layer == "l0")
                return cuts.IsInside("l0", silE, Qave);
            else if(layer == "r0")
                return cuts.IsInside("r0", silE, Qave);
        },
        {"Layer", "SilESelectedParticle", "QaveSelectedParticle"});

    TCanvas* c0 = new TCanvas("Pipe1_0", "Index of particle hitting lateral silicon", 800, 600);
    hIdx->DrawClone();

    TCanvas* c1 = new TCanvas("Pipe1_1", "PID plot", 1600, 1200);
    c1->Divide(3, 1);
    c1->cd(1);
    hPIDl0->DrawClone("colz");
    cuts.DrawCut("l0");
    c1->cd(2);
    hPIDr0->DrawClone("colz");
    cuts.DrawCut("r0");
    c1->cd(3);
    hPIDf0->DrawClone("colz");
    cuts.DrawCut("f0");

    auto dfOverlap = dfLayer.Define("SilLayerMask",
                                    [&](ActRoot::SilData& silIn)
                                    {
                                        ActRoot::SilData sil {silIn};
                                        sil.ApplyFinerThresholds(sils);

                                        bool hasF0 = sil.fSiN.count("f0") && !sil.fSiN.at("f0").empty();
                                        bool hasL0 = sil.fSiN.count("l0") && !sil.fSiN.at("l0").empty();
                                        bool hasR0 = sil.fSiN.count("r0") && !sil.fSiN.at("r0").empty();

                                        // bitmask: bit0=f0, bit1=l0, bit2=r0
                                        int mask = (hasF0 ? 1 : 0) | (hasL0 ? 2 : 0) | (hasR0 ? 4 : 0);
                                        return mask;
                                    },
                                    {"SilData"});
    std::ofstream out1(TString::Format("./Outputs/Pipe1_PIDM4_F0andSide_butNotSideLayerCount.dat").Data());
    dfOverlap.Foreach(
        [&](ActRoot::MergerData& m, const std::string& layer, int mask)
        {
            if(layer == "f0" && (mask == 3 || mask == 5)) // f0 and l0 or r0 but not both
                m.Stream(out1);
        },
        {"MergerData", "Layer", "SilLayerMask"});
    out1.close();

    auto hOverlap =
        dfOverlap.Histo1D({"hOverlap", "Silicon layer co-hits;Category;Counts", 8, -0.5, 7.5}, "SilLayerMask");

    TCanvas* cOverlap = new TCanvas("cOverlap", "f0 vs l0/r0 overlap", 800, 600);
    hOverlap->DrawClone("hist text0");

    // Etiquetas legibles para cada bin
    const char* labels[8] = {"None", "F0", "L0", "F0+L0", "R0", "F0+R0", "L0+R0", "F0+L0+R0"};
    for(int i = 0; i < 8; ++i)
        hOverlap->GetXaxis()->SetBinLabel(i + 1, labels[i]);
    hOverlap->DrawClone("hist text0");

    std::cout << "Counts per category:\n";
    for(int i = 1; i <= 8; ++i)
        std::cout << "  " << hOverlap->GetXaxis()->GetBinLabel(i) << ": " << hOverlap->GetBinContent(i) << '\n';

    std::cout << "Raw con hit lateral (l0/r0, con o sin f0): 643 (ya calculado)\n";
    std::cout
        << "Tras geometric matching (Layer != \"\"): "
        << dfWithIndex.Filter([](const std::string& l) { return l == "l0" || l == "r0"; }, {"Layer"}).Count().GetValue()
        << std::endl;
    std::cout << "Tras corte PID: " << dfPIDfiltered.Count().GetValue() << std::endl;

    std::ofstream goodLateral(TString::Format("./Outputs/Pipe1_PIDM4_goodlateral.dat").Data());
    dfL0.Foreach([&](ActRoot::MergerData& m) { m.Stream(goodLateral); }, {"MergerData"});
    dfR0.Foreach([&](ActRoot::MergerData& m) { m.Stream(goodLateral); }, {"MergerData"});

    std::vector<std::pair<int, int>> goodKeys;

    dfL0.Foreach([&](ActRoot::MergerData& m) { goodKeys.push_back({m.fRun, m.fEntry}); }, {"MergerData"});
    dfR0.Foreach([&](ActRoot::MergerData& m) { goodKeys.push_back({m.fRun, m.fEntry}); }, {"MergerData"});

    std::ofstream f0AndLateral_NotInPID(TString::Format("./Outputs/Pipe1_PIDM4_f0AndLateral_NotInPID.dat").Data());
    dfOverlap.Foreach(
        [&](ActRoot::MergerData& m, int mask)
        {
            if(mask == 2 || mask == 3 || mask == 4 || mask == 5) // l0 or r0 or f0+l0 or f0+r0
            {

                bool found = false;
                int i = 0;
                while(i < (int)goodKeys.size())
                {
                    if(goodKeys[i].first == m.fRun && goodKeys[i].second == m.fEntry)
                    {
                        found = true;
                        break;
                    }
                    i++;
                }

                if(!found)
                {
                    m.Stream(f0AndLateral_NotInPID);
                }
            }
        },
        {"MergerData", "SilLayerMask"});
    f0AndLateral_NotInPID.close();

    std::ofstream out(TString::Format("./Outputs/Pipe1_PIDM4_side_%s.dat", light.c_str()).Data());
    dfWithQave.Foreach(
        [&](ActRoot::MergerData& m, const std::string& layer)
        {
            if(layer != "f0")
                m.Stream(out);
        },
        {"MergerData", "Layer"});
    out.close();

    // Save dataframe in a .root file
    TString outfile = TString::Format("./Outputs/PIDM4_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str());
    dfPIDfiltered.Snapshot("Final_Tree", outfile);
    std::cout << "Saving Final_Tree in " << outfile << " from Pipe1_PIDM4" << '\n';
}
#endif