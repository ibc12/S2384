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

    // Get the L1 events
    auto lambdaIsL1 {[](ActRoot::MergerData& mer, ActRoot::ModularData& mod)
                     { return (mod.Get("GATCONF") == 8) && (mer.fLightIdx != -1); }};
    auto dfL1 = df.Filter(lambdaIsL1, {"MergerData", "ModularData"});

    // Filter events to avoid bad events for heavy or L1 analysis
    auto dfFilter = dfL1.Filter( // Check for heavier clusters than Li
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

    // Filter low charge deposition (high energy light particles, RANSAC needed)
    auto dfFilterL1 = dfFilter.Filter(
        [](ActRoot::MergerData& mer)
        {
            if(mer.fLight.fQave > 300.0)
                return true;
            else
                return false;
        },
        {"MergerData"});

    int nBinsS = 20;
    double sMin = -10;
    double sMax = 10;
    int minVoxelsPerSlice = 5;
    double ds = (sMax - sMin) / nBinsS;

    auto dfSigmaLight = dfFilterL1.Define("sigmaTransS",
                                          [&](ActRoot::TPCData& tpc, ActRoot::MergerData& mer)
                                          {
                                              std::vector<std::pair<double, double>> out;

                                              int idx = mer.fLightIdx;
                                              if(idx < 0)
                                                  return out;

                                              auto cluster = tpc.fClusters.at(idx);
                                              auto line = cluster.GetLine();
                                              line.Scale(2, driftFactor); // Scale manually, because we can't use function ScalVoxels

                                              auto u = line.GetDirection().Unit();
                                              auto p0 = line.GetPoint(); // centroide de carga

                                              for(int ib = 0; ib < nBinsS; ++ib)
                                              {
                                                  double sLow = sMin + ib * ds;
                                                  double sHigh = sLow + ds;
                                                  double sCenter = sLow + ds / 2.0;

                                                  double sumQ = 0.0;
                                                  double sumR2 = 0.0;
                                                  int nVox = 0;

                                                  for(auto& vox : cluster.GetVoxels())
                                                  {
                                                      for(auto& vext : vox.GetExtended())
                                                      {
                                                          vext.SetPosition({vext.GetPosition().X() * 2, vext.GetPosition().Y() * 2, vext.GetPosition().Z() * driftFactor});
                                                          auto dr = vext.GetPosition() - mer.fRP;
                                                          double s = dr.Dot(u);

                                                          if(s < sLow || s >= sHigh)
                                                              continue;

                                                          double q = vext.GetCharge();

                                                          auto d = vext.GetPosition() - p0;
                                                          double par = d.Dot(u);
                                                          double r2 = d.Mag2() - par * par;

                                                          sumQ += q;
                                                          sumR2 += q * r2;
                                                          nVox++;
                                                      }
                                                  }

                                                  // min ammount of voxels and charge to consider the slice
                                                  if(nVox < minVoxelsPerSlice || sumQ <= 0.)
                                                      continue;

                                                  double sigma = std::sqrt(sumR2 / sumQ);
                                                  out.emplace_back(sCenter, sigma);
                                              }

                                              return out;
                                          },
                                          {"TPCData", "MergerData"});


    auto hSigmaS = new TH2D("hSigmaS", "#sigma_{trans} vs s (light);s from RP [mm];#sigma_{trans} [mm]", nBinsS, sMin,
                            sMax, 100, 0, 5);

    dfSigmaLight.Foreach(
        [&](const std::vector<std::pair<double, double>>& v)
        {
            for(const auto& [s, sigma] : v)
                hSigmaS->Fill(s, sigma);
        },
        {"sigmaTransS"});

    auto prof = hSigmaS->ProfileX();


    // auto hSigmaTheta = dfSigma.Histo2D(
    //     {"hTheta", "Sigma vs ThetaHeavy;#sigma_{trans} [mm];#theta_{Heavy} [deg]", 500, 0, 5, 100, 0, 50},
    //     "sigmaTrans", "fThetaHeavy");
    // auto hSigmaPhi =
    //     dfSigma.Histo2D({"hPhi", "Sigma vs PhiHeavy;#sigma_{trans} [mm];#phi_{Heavy} [deg]", 50, 0, 4, 360, -180,
    //     180},
    //                     "sigmaTrans", "fPhiHeavy");

    // Get the histograms to avooid RDataFrame smart pointer issues
    // auto* hTheta = hSigmaTheta.GetPtr();
    // auto* hPhi = hSigmaPhi.GetPtr();

    ///////////
    // CANVAS
    ///////////

    // auto* cSigmaDebug = new TCanvas("cSigmaTheta", "Sigma vs Theta", 800, 600);
    // cSigmaDebug->Divide(2, 1);
    // cSigmaDebug->cd(1);
    // hPhi->DrawClone("COLZ");
    // cSigmaDebug->cd(2);
    // hTheta->DrawClone("COLZ");

    auto* c1 = new TCanvas("c1", "Transverse sigma of heavy particle", 800, 600);
    c1->Divide(2, 1);
    c1->cd(1);
    hSigmaS->DrawClone();
    c1->cd(2);
    prof->DrawClone();
    // Write mean and sigma in canvas
    // TLatex latex;
    // latex.SetNDC();
    // latex.SetTextSize(0.03);
    // latex.DrawLatex(0.6, 0.8, TString::Format("Mean = %.2f mm", mean));
    // latex.DrawLatex(0.6, 0.75, TString::Format("#sigma = %.2f mm", sigma));

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
}