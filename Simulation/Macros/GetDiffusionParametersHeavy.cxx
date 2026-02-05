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
#include "TF1.h"
#include "TH2.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TString.h"

#include <fstream>
#include <iostream>
#include <map>
#include <string>

#include "./../../PrettyStyle.C"

void GetDiffusionParametersHeavy()
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
    auto dfHeavy =
        dfHeavyFilter.Define("sigmaTrans",
                             [&](ActRoot::TPCData& tpc, ActRoot::MergerData& mer)
                             {
                                 int idx = mer.fHeavyIdx;
                                 auto cluster = tpc.fClusters.at(idx); // copia segura
                                 // cluster.ScaleVoxels(2, driftFactor);

                                 // --- total charge ---
                                 double sumQ = 0.;

                                 for(const auto& vox : cluster.GetVoxels())
                                 {
                                     double q = vox.GetCharge();
                                     sumQ += q;
                                 }
                                 if(sumQ == 0.)
                                     return -1.0;

                                 // --- direction ---
                                 auto line = cluster.GetLine();
                                 line.Scale(2, driftFactor); // Scale manually, because we can't use function ScalVoxels
                                 auto u = line.GetDirection().Unit();
                                 auto point = line.GetPoint(); // Charge wieghted centroid, always part of line

                                 // --- <r_trans^2> ---
                                 double sumR2 = 0.;

                                 for(auto& vox : cluster.GetVoxels())
                                 {
                                     for(auto& vext : vox.GetExtended())
                                     {
                                         vext.SetPosition({vext.GetPosition().X() * 2, vext.GetPosition().Y() * 2,
                                                           vext.GetPosition().Z() * driftFactor});
                                         double q = vext.GetCharge();
                                         auto d = vext.GetPosition() - point;
                                         double par = d.Dot(u);
                                         double r2 = d.Mag2() - par * par;
                                         sumR2 += q * r2;
                                     }
                                 }

                                 return std::sqrt(sumR2 / sumQ);
                             },
                             {"TPCData", "MergerData"});


    // Check if GetPosition gives center of voxel or corner
    auto dfVoxel = dfHeavy.Define("voxelPos",
                                  [](ActRoot::TPCData& tpc, ActRoot::MergerData& mer)
                                  {
                                      int idx = mer.fHeavyIdx;
                                      auto cluster = tpc.fClusters.at(idx); // copia segura
                                      auto voxels = cluster.GetVoxels();
                                      if(voxels.size() > 0)
                                          return voxels[0].GetPosition().Y();
                                      else
                                          return -9999.0f;
                                  },
                                  {"TPCData", "MergerData"});

    auto hVoxelPos = dfVoxel.Histo1D({"hVoxelPos", "Voxel position Y;Y [mm];Counts", 800, 0, 200}, "voxelPos");
    auto* hVoxelPosPtr = hVoxelPos.GetPtr();

    auto hSigmaTheta = dfHeavy.Histo2D(
        {"hTheta", "Sigma vs ThetaHeavy;#sigma_{trans} [mm];#theta_{Heavy} [deg]", 500, 0, 5, 100, 0, 50}, "sigmaTrans",
        "fThetaHeavy");
    auto hSigmaPhi =
        dfHeavy.Histo2D({"hPhi", "Sigma vs PhiHeavy;#sigma_{trans} [mm];#phi_{Heavy} [deg]", 50, 0, 4, 360, -180, 180},
                        "sigmaTrans", "fPhiHeavy");

    // Get the histograms to avooid RDataFrame smart pointer issues
    auto* hTheta = hSigmaTheta.GetPtr();
    auto* hPhi = hSigmaPhi.GetPtr();

    // Plot sigmaTrans
    auto hSigmaTrans = dfHeavy.Histo1D(
        {"hSigmaTrans", "Transverse sigma of heavy particle;#sigma_{trans} [mm];Counts", 100, 0, 3}, "sigmaTrans");
    auto* hSigmaTransPtr = hSigmaTrans.GetPtr();
    // Fit to gausian in the range 0-7 mm
    hSigmaTransPtr->Fit("gaus", "", "", 0, 3);
    // Get sigma and mean of fit to put in canvas
    auto fit = hSigmaTransPtr->GetFunction("gaus");
    double mean = fit->GetParameter(1);
    double sigma = fit->GetParameter(2);

    // Plot voxel position histogram
    TCanvas* cVoxel = new TCanvas("cVoxel", "Voxel Position Y", 800, 600);
    cVoxel->cd();
    hVoxelPosPtr->DrawClone();


    ///////////
    // CANVAS
    ///////////

    auto* cSigmaDebug = new TCanvas("cSigmaTheta", "Sigma vs Theta", 800, 600);
    cSigmaDebug->Divide(2, 1);
    cSigmaDebug->cd(1);
    hPhi->DrawClone("COLZ");
    cSigmaDebug->cd(2);
    hTheta->DrawClone("COLZ");

    auto* c1 = new TCanvas("c1", "Transverse sigma of heavy particle", 800, 600);
    c1->cd();
    hSigmaTransPtr->DrawClone();
    // Write mean and sigma in canvas
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.03);
    latex.DrawLatex(0.6, 0.8, TString::Format("Mean = %.2f mm", mean));
    latex.DrawLatex(0.6, 0.75, TString::Format("#sigma = %.2f mm", sigma));

    // Look at the events to ensure they are good candidates
    // std::ofstream outFile {"./Outputs/Events_High_SigmaT.dat"};
    // dfHeavy.Foreach(
    //     [&outFile](const ActRoot::MergerData& mer, double sigmaT)
    //     {
    //         if(sigmaT > 1.7)
    //             mer.Stream(outFile);
    //     },
    //     {"MergerData", "sigmaTrans"});
    // outFile.close();
    //
    // std::ofstream outFile1 {"./Outputs/Events_Low_SigmaT.dat"};
    // dfHeavy.Foreach(
    //     [&outFile1](const ActRoot::MergerData& mer, double sigmaT)
    //     {
    //         if(sigmaT <= 1.2)
    //             mer.Stream(outFile1);
    //     },
    //     {"MergerData", "sigmaTrans"});
    // outFile1.close();
}