#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActInputParser.h"
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

#include <fstream>
#include <iostream>
#include <map>
#include <string>

#include "./../../PrettyStyle.C"

void GetDiffusionParameters()
{
    // Get diffusion parameters from simulation
    // PrettyStyle(false);
    std::string dataconf {"../../configs/data.conf"};

    // Read config file to get drift parameter
    ActRoot::InputParser parser {};
    parser.ReadFile("../../configs/detector.conf");
    auto driftBlock = parser.GetBlock("Merger");
    auto driftFactor = driftBlock->GetDouble("DriftFactor"); // in mm^2/us


    // Read data
    ActRoot::DataManager dataman {dataconf, ActRoot::ModeType::EMerge};
    dataman.SetRuns(64, 67);
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    auto chain3 {dataman.GetChain(ActRoot::ModeType::EFilter)};
    auto chain4 {dataman.GetChain(ActRoot::ModeType::EReadTPC)};
    chain->AddFriend(chain2.get());
    chain->AddFriend(chain3.get(), "TPCData");
    chain->AddFriend(chain4.get(), "GETTree");

    // RDataFrame
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {*chain}; 

    // df.Describe().Print();

    // Filter events to avoid bad events for heavy analysis
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
                            {"MergerData", "GETTree.TPCData"});

    // Get the L1 events
    auto lambdaIsL1 {[](ActRoot::MergerData& mer, ActRoot::ModularData& mod)
                     { return (mod.Get("GATCONF") == 8) && (mer.fLightIdx != -1); }};
    auto dfL1 = dfFilter.Filter(lambdaIsL1, {"MergerData", "ModularData"});

    // Go for events with z almost continuous
    // Filter low charge deposition (high energy light particles, RANSAC needed)
    auto dfFilterL1 = dfL1.Filter(
                              [](ActRoot::MergerData& mer, ActRoot::TPCData& tpc)
                              {
                                  auto thetaLight = mer.fThetaLight;
                                  auto phiLight = std::abs(mer.fPhiLight); // aboid -90 deg ambiguity
                                  if(phiLight < 100.0 && phiLight > 80.0 && mer.fLightIdx != -1)
                                      return true;
                                  else if(phiLight > 260.0 && phiLight < 280.0 && mer.fLightIdx != -1)
                                      return true;
                                  else
                                      return false;
                              },
                              {"MergerData", "TPCData"})
                          .Filter(
                              [](ActRoot::MergerData& mer)
                              {
                                  if(mer.fLight.fQave > 300.0)
                                      return true;
                                  else
                                      return false;
                              },
                              {"MergerData"});

    // Debug with heavy particle
    auto dfHeavyFilter = dfFilter.Filter(
        [](ActRoot::MergerData& mer)
        {
            if(mer.fHeavyIdx != -1)
                return true;
            else
                return false;
        },
        {"MergerData"});
    auto dfHeavy = dfHeavyFilter.Define("sigmaTrans",
                                        [&](ActRoot::TPCData& tpc, ActRoot::MergerData& mer)
                                        {
                                            int idx = mer.fHeavyIdx;
                                            auto cluster = tpc.fClusters.at(idx); // copia segura
                                            cluster.ScaleVoxels(2, driftFactor);

                                            // --- centro de carga ---
                                            ROOT::Math::XYZVector rmean(0., 0., 0.);
                                            double sumQ = 0.;

                                            for(const auto& vox : cluster.GetVoxels())
                                            {
                                                double q = vox.GetCharge();
                                                rmean += q * ROOT::Math::XYZVector(vox.GetPosition());
                                                sumQ += q;
                                            }

                                            if(sumQ == 0.)
                                                return -1.0;

                                            rmean *= (1. / sumQ);

                                            // --- dirección ---
                                            auto u = cluster.GetLine().GetDirection().Unit();

                                            // --- <r_perp^2> ---
                                            double sumR2 = 0.;

                                            for(const auto& vox : cluster.GetVoxels())
                                            {
                                                double q = vox.GetCharge();
                                                auto d = vox.GetPosition() - rmean;
                                                double par = d.Dot(u);
                                                double r2 = d.Mag2() - par * par;
                                                sumR2 += q * r2;
                                            }

                                            return std::sqrt(sumR2 / sumQ);
                                        },
                                        {"TPCData", "MergerData"});


    auto hSigmaTheta =
        dfHeavy.Histo2D({"h", "Sigma vs ThetaHeavy", 500, 0, 10, 100, 0, 50}, "sigmaTrans", "fThetaHeavy");
    auto hSigmaPhi = dfHeavy.Histo2D({"h", "Sigma vs PhiHeavy", 50, 0, 4, 360, -180, 180}, "sigmaTrans", "fPhiHeavy");
    auto* cSigmaDebug = new TCanvas("cSigmaTheta", "Sigma vs Theta", 800, 600);
    cSigmaDebug->Divide(2, 1);
    cSigmaDebug->cd(1);
    hSigmaPhi->DrawClone("COLZ");
    cSigmaDebug->cd(2);
    hSigmaTheta->DrawClone("COLZ");

    // Debug unit test define and histogram last voxel of heavy particle
    auto dfTest = dfHeavy.Define("lastVoxelX",
                                 [&](ActRoot::TPCData& tpc, ActRoot::MergerData& mer)
                                 {
                                     int heavyIdx = mer.fHeavyIdx;
                                     auto cluster = tpc.fClusters.at(heavyIdx);
                                     cluster.ScaleVoxels(2, driftFactor); // scale according to pad size and drift
                                     cluster.SortAlongDir(cluster.GetLine().GetDirection().Unit());
                                     auto voxels = cluster.GetRefToVoxels();
                                     if(voxels.size() == 0)
                                         return -1.0f;
                                     auto lastVoxel = voxels.back();
                                     float position = lastVoxel.GetPosition().X();
                                     return position;
                                 },
                                 {"TPCData", "MergerData"});

    auto hLastVoxelX = dfTest.Histo1D(
        {"hLastVoxelX", "X position of last voxel of heavy particle;X [mm];Counts", 100, 0, 300}, "lastVoxelX");
    auto* c0 = new TCanvas("c0", "X position of last voxel of heavy particle", 800, 600);
    c0->cd();
    hLastVoxelX->DrawClone();

    // Debug line point position for heavy particle
    auto dfTestLinePoint =
        dfHeavy.Define("linePointX",
                       [&](ActRoot::MergerData& mer, ActRoot::TPCData& tpc)
                       {
                           int heavyIdx = mer.fHeavyIdx;
                           auto cluster = tpc.fClusters.at(heavyIdx);
                           cluster.ScaleVoxels(2, driftFactor); // scale according to pad size and drift
                           auto r0 = cluster.GetRefToLine().GetPoint();
                           return r0.X();
                       },
                       {"MergerData", "TPCData"});
    auto hLinePointX = dfTestLinePoint.Histo1D(
        {"hLinePointX", "X position of line point of heavy particle;X [mm];Counts", 100, -100, 400}, "linePointX");
    auto* cLinePoint = new TCanvas("cLinePoint", "X position of line point of heavy particle", 800, 600);
    cLinePoint->cd();
    hLinePointX->DrawClone();


    // Plot sigmaTrans
    auto hSigmaTrans = dfHeavy.Histo1D(
        {"hSigmaTrans", "Transverse sigma of heavy particle;#sigma_{trans} [mm];Counts", 100, 0, 3}, "sigmaTrans");
    auto* c1 = new TCanvas("c1", "Transverse sigma of heavy particle", 800, 600);
    c1->cd();
    hSigmaTrans->DrawClone();

    // Look at the events to ensure they are good candidates
    // std::ofstream outFile {"./Outputs/Events_High_SigmaT.dat"};
    // dfHeavy.Foreach(
    //     [&outFile](const ActRoot::MergerData& mer, double sigmaT)
    //     {
    //         if(sigmaT > 7)
    //             mer.Stream(outFile);
    //     },
    //     {"MergerData", "sigmaTrans"});
    // outFile.close();
    // 
    // std::ofstream outFile1 {"./Outputs/Events_Low_SigmaT.dat"};
    // dfHeavy.Foreach(
    //     [&outFile1](const ActRoot::MergerData& mer, double sigmaT)
    //     {
    //         if(sigmaT <= 2)
    //             mer.Stream(outFile1);
    //     },
    //     {"MergerData", "sigmaTrans"});
    // outFile1.close();
}