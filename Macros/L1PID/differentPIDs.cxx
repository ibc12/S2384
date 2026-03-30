#include "ActCluster.h"
#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActLine.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActSRIM.h"
#include "ActTPCData.h"
#include "ActTPCParameters.h"
#include "ActVoxel.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"
#include <random>

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMath.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TSpline.h>
#include <TStyle.h>

#include <Math/Point3D.h>
#include <Math/Vector3D.h>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

#include "../../PrettyStyle.C"


void differentPIDs()
{
    PrettyStyle(true, true);

    // Get all L1 experimental data for 7Li; cut in L1 and filter bad events, select elastic, check good fits to
    // deuterium

    // Get all data for 7Li (there is no triton there, so easier to see deuterium)
    std::string beam {"7Li"};
    std::string target {"d"};
    std::string light {"p"};

    std::string dataconf {};
    if(beam == "11Li")
        dataconf = "../../configs/data_11Li.conf";
    else if(beam == "7Li")
        dataconf = "../../configs/data_7Li.conf";
    else
        throw std::runtime_error("Unknown beam");

    // Read data
    ActRoot::DataManager dataman {dataconf, ActRoot::ModeType::EMerge};
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    auto chain3 {dataman.GetChain(ActRoot::ModeType::EFilter)};
    chain->AddFriend(chain2.get());
    chain->AddFriend(chain3.get());

    // Get drift parameter
    ActRoot::InputParser parser {};
    parser.ReadFile("../../configs/detector.conf");
    auto driftBlock = parser.GetBlock("Merger");
    auto driftFactor = driftBlock->GetDouble("DriftFactor"); // in mm^2/us

    // RDataFrame
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame dforigin {*chain};

    // Filter GATCONF == L1
    auto df =
        dforigin.Filter([](ActRoot::ModularData& mod, ActRoot::MergerData& mer)
                        { return mod.Get("GATCONF") == 8 && (mer.fLightIdx != -1); }, {"ModularData", "MergerData"});

    auto dfFilter = df.Filter( // Check for heavier clusters than Li
                          [](ActRoot::MergerData& m)
                          {
                              if(m.fHeavy.fQave > 2000.)
                                  return false;
                              return true;
                          },
                          {"MergerData"})
                        .Filter( // Check if heavy reaches end of ACTAR (other way to mask heavier clusters,
                                 // but can delete good events with heavy particle bad reconstructed)
                            [](ActRoot::MergerData& m)
                            {
                                auto TL {m.fHeavy.fTL};
                                auto theta {m.fThetaHeavy * TMath::DegToRad()};
                                auto z_end {m.fRP.X() + TL * TMath::Cos(theta)};
                                if(z_end < 240)
                                    return false;
                                return true;
                            },
                            {"MergerData"})
                        .Filter( // Most of the times high charge deposit
                                 // is masked by rp or is in beam cluster
                            [&](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                            {
                                auto rp {m.fRP};
                                auto rp_y {rp.Y() / 2};
                                // Run for all clusters
                                int counter {};
                                for(auto& cluster : tpc.fClusters)
                                {
                                    auto voxels {cluster.GetRefToVoxels()};
                                    for(auto& v : voxels)
                                    {
                                        if(v.GetPosition().Y() > rp_y - 3 &&
                                           v.GetPosition().Y() < rp_y + 3) // aprox L1 exclusion zone
                                            if(v.GetCharge() > 3000.)
                                                counter++;
                                    }
                                }
                                if(counter > 4)
                                    return false;
                                return true;
                            },
                            {"MergerData", "TPCData"});

    auto dfZandQ = dfFilter.Define("zDrift",
                                   [&](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                                   {
                                       auto& rp {m.fRP};

                                       auto idx = m.fLightIdx;
                                       if(idx < 0)
                                           return -100.0f;

                                       if(static_cast<std::size_t>(idx) >= tpc.fClusters.size())
                                           return -100.0f;

                                       auto& cluster = tpc.fClusters.at(idx);
                                       auto& voxels = cluster.GetRefToVoxels();

                                       if(tpc.fRPs.empty())
                                           return -100.0f;
                                       auto rpVox = tpc.fRPs.front();

                                       float maxDist = -1.0;
                                       float zExtreme = 0.0;

                                       for(auto& v : voxels)
                                       {
                                           auto pos = v.GetPosition();

                                           float dx = pos.X() - rpVox.X();
                                           float dy = pos.Y() - rpVox.Y();
                                           float dz = pos.Z() - rpVox.Z();

                                           float dist = std::sqrt(dx * dx + dy * dy + dz * dz);

                                           if(dist > maxDist)
                                           {
                                               maxDist = dist;
                                               zExtreme = pos.Z();
                                           }
                                       }

                                       // diferencia en Z respecto al RP voxel
                                       float deltaZ = (zExtreme - rpVox.Z()) * driftFactor;

                                       float zDrift = 110.0 + deltaZ;

                                       return zDrift;
                                   },
                                   {"MergerData", "TPCData"});

    // Filter the z distances with file zDrift_cut_perRun.dat
    // Read per-run zDrift cuts: columns are run, meanStart, sigmaStart, meanEnd, sigmaEnd
    struct ZDriftCut
    {
        double zMin; // meanStart + nSigma * sigmaStart
        double zMax; // meanEnd   - nSigma * sigmaEnd
    };
    const double nSigmaZDrift {3.0};
    std::map<int, ZDriftCut> zDriftCuts;
    {
        std::ifstream finZ("./Outputs/zDrift_cut_perRun.dat");
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

    auto dfFilterZandQ = dfZandQ.Filter(
        [&zDriftCuts](float zDrift, ActRoot::MergerData& m)
        {
            auto it = zDriftCuts.find(m.fRun);
            if(it == zDriftCuts.end())
                return false; // run not in file → reject
            return zDrift > it->second.zMin && zDrift < it->second.zMax && m.fLight.fQave > 600;
        },
        {"zDrift", "MergerData"});

    auto dfVoxelLight = dfFilterZandQ.Define("chargePerVoxelLight",
                                             [&](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                                             {
                                                 auto idx = m.fLightIdx;
                                                 if(idx < 0)
                                                     return -1.0f;
                                                 auto& voxels = tpc.fClusters[idx].GetRefToVoxels();
                                                 float totalCharge = 0.0f;
                                                 for(auto& v : voxels)
                                                     totalCharge += v.GetCharge();
                                                 return totalCharge / voxels.size();
                                             },
                                             {"MergerData", "TPCData"});

    // Get the cut for PID
    ActRoot::CutsManager<std::string> cuts;
    cuts.ReadCut("l1", TString::Format("./Cuts/%s_TLvsQ_%s.root", light.c_str(), beam.c_str()).Data());
    cuts.ReadCut("l1_theta", TString::Format("./Cuts/%s_ThetaVSq_%s.root", light.c_str(), beam.c_str()).Data());
    cuts.ReadCut("l1_chargePerVoxel",
                 TString::Format("./Cuts/%s_TLvsQvoxels_%s.root", light.c_str(), beam.c_str()).Data());

    // Apply cut in theta - Qtotal
    auto dfLight = dfVoxelLight.Filter(
        [&](ActRoot::MergerData& m) { return cuts.IsInside("l1", m.fLight.fRawTL, m.fLight.fQtotal); }, {"MergerData"});
    auto dfLightAngle =
        dfLight.Filter([&](ActRoot::MergerData& m)
                       { return cuts.IsInside("l1_theta", m.fThetaLight, m.fLight.fQtotal); }, {"MergerData"});

    auto dfVoxelLightCut =
        dfVoxelLight.Filter([&](ActRoot::MergerData& m, float chargePerVoxelLight)
                            { return cuts.IsInside("l1_chargePerVoxel", m.fLight.fRawTL, chargePerVoxelLight); },
                            {"MergerData", "chargePerVoxelLight"});

    auto dfVoxelLightCut_PIDcut = dfVoxelLightCut.Filter(
        [&](ActRoot::MergerData& m) { return cuts.IsInside("l1", m.fLight.fRawTL, m.fLight.fQtotal); }, {"MergerData"});


    // PIDs for filter Z and Q data
    auto hTL_Qtot_filterZandQ = dfVoxelLight.Histo2D(
        {"hTL_Qtot_filterZandQ", "Track Length vs Qtotal with Z and Q cut;Track Length [a.u.];Qtotal [a.u.]", 240, 0,
         120, 2000, 0, 3e5},
        "MergerData.fLight.fRawTL", "MergerData.fLight.fQtotal");
    auto hTheta_Qtot_filterZandQ =
        dfVoxelLight.Histo2D({"hTheta_Qtot_filterZandQ", "Theta vs Qtotal with Z and Q cut;Theta [deg];Qtotal [a.u.]",
                              180, 0, 180, 2000, 0, 3e5},
                             "MergerData.fThetaLight", "MergerData.fLight.fQtotal");
    // Now after the cut in TL-Q
    auto hTL_Qtot_cut =
        dfLight.Histo2D({"hTL_Qtot_cut", "Track Length vs Qtotal with cut;Track Length [a.u.];Qtotal [a.u.]", 240, 0,
                         120, 2000, 0, 3e5},
                        "MergerData.fLight.fRawTL", "MergerData.fLight.fQtotal");
    auto hTheta_Qtot_cut = dfLight.Histo2D(
        {"hTheta_Qtot_cut", "Theta vs Qtotal with cut;Theta [deg];Qtotal [a.u.]", 180, 0, 180, 2000, 0, 3e5},
        "MergerData.fThetaLight", "MergerData.fLight.fQtotal");

    // PID plots per voxel charge
    auto hTL_Qtot_chargePerVoxel = dfVoxelLight.Histo2D(
        {"hTL_Qtot_chargePerVoxel", "Track Length vs Qtotal per voxel;Track Length [a.u.];Qtotal [a.u.]", 240, 0, 120,
         2000, 0, 3e3},
        "MergerData.fLight.fRawTL", "chargePerVoxelLight");
    auto hTheta_Qtot_chargePerVoxel =
        dfVoxelLight.Histo2D({"hTheta_Qtot_chargePerVoxel", "Theta vs Qtotal per voxel;Theta [deg];Qtotal [a.u.]", 180,
                              0, 180, 2000, 0, 3e3},
                             "MergerData.fThetaLight", "chargePerVoxelLight");
    auto hTL_Qtot_cut_chargePerVoxel =
        dfLight.Histo2D({"hTL_Qtot_cut_chargePerVoxel",
                         "Track Length vs Qtotal with cut in original PID;Track Length [a.u.];Qtotal [a.u.]", 240, 0,
                         120, 2000, 0, 3e3},
                        "MergerData.fLight.fRawTL", "chargePerVoxelLight");
    auto hTheta_Qtot_cut_chargePerVoxel = dfLight.Histo2D(
        {"hTheta_Qtot_cut_chargePerVoxel", "Theta vs Qtotal with cut in original PID;Theta [deg];Qtotal [a.u.]", 180, 0,
         180, 2000, 0, 3e3},
        "MergerData.fThetaLight", "chargePerVoxelLight");

    // PID plots for cut in TL-Q per voxel
    auto hTL_Qtot_cut_chargePerVoxel_cut_CUT = dfVoxelLightCut.Histo2D(
        {"hTL_Qtot_cut_chargePerVoxel_cut_CUT",
         "Track Length vs Qtotal per voxel with cut by charge per voxel;Track Length [a.u.];Qtotal [a.u.]", 240, 0, 120,
         2000, 0, 3e3},
        "MergerData.fLight.fRawTL", "chargePerVoxelLight");
    auto hTheta_Qtot_cut_chargePerVoxel_cut_CUT = dfVoxelLightCut.Histo2D(
        {"hTheta_Qtot_cut_chargePerVoxel_cut_CUT",
         "Theta vs Qtotal per voxel with cut by charge per voxel;Theta [deg];Qtotal [a.u.]", 180, 0, 180, 2000, 0, 3e3},
        "MergerData.fThetaLight", "chargePerVoxelLight");
    auto hTL_Qtot_cut_chargePerVoxel_CUT = dfVoxelLightCut.Histo2D(
        {"hTL_Qtot_cut_chargePerVoxel_CUT",
         "Track Length vs Qtotal with cut by charge per voxel;Track Length [a.u.];Qtotal [a.u.]", 240, 0, 120, 2000, 0,
         3e5},
        "MergerData.fLight.fRawTL", "MergerData.fLight.fQtotal");
    auto hTheta_Qtot_cut_chargePerVoxel_CUT = dfVoxelLightCut.Histo2D(
        {"hTheta_Qtot_cut_chargePerVoxel_CUT", "Theta vs Qtotal with cut by charge per voxel;Theta [deg];Qtotal [a.u.]",
         180, 0, 180, 2000, 0, 3e5},
        "MergerData.fThetaLight", "MergerData.fLight.fQtotal");

    // PID plots for cut in TL-Q per voxel and original PID cut
    auto hTL_Qtot_cut_chargePerVoxel_cut_CUT_PIDcut =
        dfVoxelLightCut_PIDcut.Histo2D({"hTL_Qtot_cut_chargePerVoxel_cut_CUT_PIDcut",
                                        "Track Length vs Qtotal per voxel with cut by charge per voxel and original "
                                        "PID cut;Track Length [a.u.];Qtotal [a.u.]",
                                        240, 0, 120, 2000, 0, 3e5},
                                       "MergerData.fLight.fRawTL", "MergerData.fLight.fQtotal");
    auto hTheta_Qtot_cut_chargePerVoxel_cut_CUT_PIDcut = dfVoxelLightCut_PIDcut.Histo2D(
        {"hTheta_Qtot_cut_chargePerVoxel_cut_CUT_PIDcut",
         "Theta vs Qtotal with cut by charge per voxel and original PID cut;Theta [deg];Qtotal [a.u.]", 180, 0, 180,
         2000, 0, 3e5},
        "MergerData.fThetaLight", "MergerData.fLight.fQtotal");
    auto hTL_Qtot_chargePerVoxel_cut_chargePerVoxel_CUT_PIDcut = dfVoxelLightCut_PIDcut.Histo2D(
        {"hTL_Qtot_cut_chargePerVoxel_CUT_PIDcut",
         "Track Length vs Qtotal with cut by charge per voxel and original PID cut;Track Length [a.u.];Qtotal [a.u.]",
         240, 0, 120, 2000, 0, 3e3},
        "MergerData.fLight.fRawTL", "chargePerVoxelLight");
    auto hTheta_Qtot_chargePerVoxel_cut_chargePerVoxel_CUT_PIDcut = dfVoxelLightCut_PIDcut.Histo2D(
        {"hTheta_Qtot_cut_chargePerVoxel_CUT_PIDcut",
         "Theta vs Qtotal with cut by charge per voxel and original PID cut;Theta [deg];Qtotal [a.u.]", 180, 0, 180,
         2000, 0, 3e3},
        "MergerData.fThetaLight", "chargePerVoxelLight");


    // Plot
    auto* c1 = new TCanvas("c1", "c1", 1200, 800);
    c1->Divide(2, 2);
    c1->cd(1);
    hTL_Qtot_filterZandQ->DrawClone("COLZ");
    cuts.DrawCut("l1");
    c1->cd(2);
    hTL_Qtot_cut->DrawClone();
    c1->cd(3);
    hTheta_Qtot_filterZandQ->DrawClone("COLZ");
    cuts.DrawCut("l1_theta");
    c1->cd(4);
    hTheta_Qtot_cut->DrawClone("COLZ");

    auto* c2 = new TCanvas("c2", "c2", 1200, 600);
    c2->Divide(2, 2);
    c2->cd(1);
    hTL_Qtot_chargePerVoxel->DrawClone("COLZ");
    cuts.DrawCut("l1_chargePerVoxel");
    c2->cd(2);
    hTheta_Qtot_chargePerVoxel->DrawClone("COLZ");
    c2->cd(3);
    hTL_Qtot_cut_chargePerVoxel->DrawClone("COLZ");
    c2->cd(4);
    hTheta_Qtot_cut_chargePerVoxel->DrawClone("COLZ");

    auto* c3 = new TCanvas("c3", "c3", 1200, 600);
    c3->Divide(2, 2);
    c3->cd(1);
    hTL_Qtot_cut_chargePerVoxel_cut_CUT->DrawClone("COLZ");
    c3->cd(2);
    hTheta_Qtot_cut_chargePerVoxel_cut_CUT->DrawClone("COLZ");
    c3->cd(3);
    hTL_Qtot_cut_chargePerVoxel_CUT->DrawClone("COLZ");
    c3->cd(4);
    hTheta_Qtot_cut_chargePerVoxel_CUT->DrawClone("COLZ");

    auto * c4 = new TCanvas("c4", "c4", 1200, 600);
    c4->Divide(2, 2);
    c4->cd(1);
    hTL_Qtot_cut_chargePerVoxel_cut_CUT_PIDcut->DrawClone("COLZ");
    c4->cd(2);
    hTheta_Qtot_cut_chargePerVoxel_cut_CUT_PIDcut->DrawClone("COLZ");
    c4->cd(3);
    hTL_Qtot_chargePerVoxel_cut_chargePerVoxel_CUT_PIDcut->DrawClone("COLZ");
    c4->cd(4);
    hTheta_Qtot_chargePerVoxel_cut_chargePerVoxel_CUT_PIDcut->DrawClone("COLZ");
}