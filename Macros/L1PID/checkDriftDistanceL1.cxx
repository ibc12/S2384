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
#include <TGraphErrors.h>
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


void checkDriftDistanceL1()
{
    PrettyStyle(true, true);

    // Get all L1 experimental data for 7Li; cut in L1 and filter bad events, select elastic, check good fits to
    // deuterium

    // Get all data for 7Li (there is no triton there, so easier to see deuterium)
    std::string beam {"7Li"};
    std::string target {"d"};
    std::string light {"p"};

    std::string dataconf {"../../configs/data.conf"};
    // if(beam == "11Li")
    //     dataconf = "../../configs/data_11Li.conf";
    // else if(beam == "7Li")
    //     dataconf = "../../configs/data_7Li.conf";
    // else
    //     throw std::runtime_error("Beam cannot differ from 11Li or 7Li");

    // Read data
    ActRoot::DataManager dataman {dataconf, ActRoot::ModeType::EMerge};
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    auto chain3 {dataman.GetChain(ActRoot::ModeType::EFilter)};
    auto chain4 {dataman.GetChain(ActRoot::ModeType::EReadTPC)};
    chain->AddFriend(chain2.get());
    chain->AddFriend(chain3.get(), "TPCData");
    chain->AddFriend(chain4.get(), "GETTree");

    // Get drift parameter
    ActRoot::InputParser parser {};
    parser.ReadFile("../../configs/detector.conf");
    auto driftBlock = parser.GetBlock("Merger");
    auto driftFactor = driftBlock->GetDouble("DriftFactor"); // in mm^2/us

    // Get SRIM files to do the fits
    auto* srim = new ActPhysics::SRIM;
    srim->ReadTable("light", "../../Calibrations/SRIM/1H_900mb_CF4_95-5.txt");
    srim->ReadTable("lightD", "../../Calibrations/SRIM/2H_900mb_CF4_95-5.txt");
    srim->ReadTable("lightT", "../../Calibrations/SRIM/3H_900mb_CF4_95-5.txt");
    srim->ReadTable("light3He", "../../Calibrations/SRIM/3He_900mb_CF4_95-5.txt");
    srim->ReadTable("light4He", "../../Calibrations/SRIM/4He_900mb_CF4_95-5.txt");

    // RDataFrame
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame dforigin {*chain};

    // Cuts
    ActRoot::CutsManager<std::string> cuts;

    // Filter GATCONF == L1
    auto df =
        dforigin.Filter([](ActRoot::ModularData& mod, ActRoot::MergerData& mer)
                        { return mod.Get("GATCONF") == 8 && (mer.fLightIdx != -1); }, {"ModularData", "MergerData"});

    // df.Describe().Print();

    // Filter bad events
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

    // auto dfZdrift = dfFilter.Define("zDrift",
    //                                 [&](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
    //                                 {
    //                                     auto& rp {m.fRP}; // This equals aprox 110 mm
    //                                     // Get last point of the track
    //                                     auto idx = m.fLightIdx;
    //                                     if(idx < 0)
    //                                         return -100.0;
    //                                     // guard against out-of-range index (this was causing vector::at to throw)
    //                                     if(static_cast<std::size_t>(idx) >= tpc.fClusters.size())
    //                                         return -100.0;
    //                                     auto& cluster = tpc.fClusters.at(idx);
    //                                     auto line = cluster.GetLine();
    //                                     auto dir = line.GetDirection().Unit();
    //                                     auto rpVox = tpc.fRPs.front();
    //                                     // cluster.SortAlongDir(line.GetDirection().Unit());
    //                                     auto& voxels = cluster.GetRefToVoxels();
    //                                     // Get mannually the z point whose difference with z rp is max
    //                                     double maxAbsDeltaZ = -1;
    //                                     double deltaZextreme = 0;
    //
    //                                     for(const auto& v : voxels)
    //                                     {
    //                                         // std::cout << "Voxel Z: " << v.GetPosition().Z()  << " RP Z: " <<
    //                                         rp.Z() << " mm" << std::endl; double z = v.GetPosition().Z() *
    //                                         driftFactor; double deltaZ = z - rp.Z();
    //
    //                                         if(std::abs(deltaZ) > maxAbsDeltaZ)
    //                                         {
    //                                             maxAbsDeltaZ = std::abs(deltaZ);
    //                                             deltaZextreme = deltaZ;
    //                                         }
    //                                     }
    //
    //                                     double zDrift = 110.0 + deltaZextreme; // mm, distance from RP to pad plane
    //                                     // auto lastZ = voxels.back().GetPosition().Z() * driftFactor;
    //                                     // Z difference relative to the RP
    //                                     //  deltaZ > 0  -> track goes to larger Z (away from pad plane)
    //                                     //  deltaZ < 0  -> track goes to smaller Z (closer to pad plane)
    //                                     double deltaZ = deltaZextreme - rp.Z();
    //                                     if(m.fRun == 66 && m.fEntry == 1740 )
    //                                     {
    //                                         std::cout << "RP Z: " << rp.Z() << " mm, deltaZextreme: " <<
    //                                         deltaZextreme << " mm, zDrift: " << zDrift << " mm" << std::endl;
    //                                     }
    //
    //                                     return zDrift;
    //                                 },
    //                                 {"MergerData", "TPCData"});
    auto dfZdrift = dfFilter.Define("zDrift",
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
                          
    auto hZdrift = dfZdrift.Histo1D({"hZdrift", "Z drift distribution;Z drift [mm];Counts", 340, -40, 300}, "zDrift");
    // Fit extremes into gaussian to get mean and sigma of the distribution
    auto fitGausStart = new TF1("fitFunc", "gaus", -20, 5);
    auto fitGausEnd = new TF1("fitFunc", "gaus", 250, 270);
    hZdrift->Fit(fitGausStart, "RQ");
    hZdrift->Fit(fitGausEnd, "RQ");
    auto meanStart = fitGausStart->GetParameter(1);
    auto sigmaStart = fitGausStart->GetParameter(2);
    auto meanEnd = fitGausEnd->GetParameter(1);
    auto sigmaEnd = fitGausEnd->GetParameter(2);
    std::cout << "Start cut in Z drift: " << meanStart + 3 * sigmaStart << " mm" << std::endl;
    std::cout << "End cut in Z drift: " << meanEnd - 3 * sigmaEnd << " mm" << std::endl;

    // Histo zDrift vs theta and Qave
    auto hZdriftTheta =
        dfZdrift.Histo2D({"hZdriftTheta", "Z drift vs #theta;#theta [deg];Z drift [mm]", 90, 0, 90, 340, -40, 300},
                         "MergerData.fThetaLight", "zDrift");
    auto hZdriftQave =
        dfZdrift.Histo2D({"hZdriftQave", "Z drift vs Qave;Qave [a.u.];Z drift [mm]", 100, 0, 2000, 340, -40, 300},
                         "MergerData.fQave", "zDrift");

    auto hQave = dfZdrift.Histo1D({"hQave", "Qave distribution;Qave [a.u.];Counts", 150, 0, 3000}, "MergerData.fQave");

    // Filter in zDrift
    auto dfZdriftFiltered = dfZdrift.Filter(
        [&](float zDrift)
        {
            if(zDrift < meanStart + 3 * sigmaStart || zDrift > meanEnd - 3 * sigmaEnd)
                return false;
            return true;
        },
        {"zDrift"});

    auto hZdriftFiltered = dfZdriftFiltered.Histo1D(
        {"hZdriftFiltered", "Z drift distribution filtered;Z drift [mm];Counts", 340, -40, 300}, "zDrift");
    auto hQaveFiltered = dfZdriftFiltered.Histo1D(
        {"hQaveFiltered", "Qave distribution filtered;Qave [a.u.];Counts", 300, 0, 3000}, "MergerData.fQave");
    auto fitGausQ = new TF1("fitFunc", "gaus", 0, 120);
    hQaveFiltered->Fit(fitGausQ, "RQ");
    auto meanQ = fitGausQ->GetParameter(1);
    auto sigmaQ = fitGausQ->GetParameter(2);
    std::cout << "Start cut in Qave: " << meanQ - 3 * sigmaQ << " a.u." << std::endl;
    std::cout << "End cut in Qave: " << meanQ + 3 * sigmaQ << " a.u." << std::endl;

    // auto cChi2 = new TCanvas("cChi2", "cChi2", 800, 600);
    // hChi2NormalizedByDeuterium->DrawClone();

    auto cZdrift = new TCanvas("cZdrift", "cZdrift", 800, 600);
    cZdrift->Divide(2, 2);
    cZdrift->cd(1);
    hZdrift->DrawClone();
    fitGausStart->SetLineColor(kRed);
    fitGausStart->Draw("same");
    fitGausEnd->SetLineColor(kBlue);
    fitGausEnd->Draw("same");
    // Write the mean and sigma of the fits in the canvas
    auto latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.04);
    latex->DrawLatex(0.15, 0.85, Form("Start fit: #mu = %.1f mm, #sigma = %.1f mm", meanStart, sigmaStart));
    latex->DrawLatex(0.15, 0.80, Form("End fit: #mu = %.1f mm, #sigma = %.1f mm", meanEnd, sigmaEnd));
    cZdrift->cd(2);
    hZdriftTheta->DrawClone("COLZ");
    cZdrift->cd(3);
    hZdriftQave->DrawClone("COLZ");
    cZdrift->cd(4);
    hQave->DrawClone();

    auto cZdriftFiltered = new TCanvas("cZdriftFiltered", "cZdriftFiltered", 800, 600);
    cZdriftFiltered->Divide(2, 1);
    cZdriftFiltered->cd(1);
    hZdriftFiltered->DrawClone();
    cZdriftFiltered->cd(2);
    hQaveFiltered->DrawClone();
    fitGausQ->SetLineColor(kRed);
    fitGausQ->Draw("same");

    // ============================================================
    // Per-run fit of zDrift extremes (pad-plane and cathode peaks)
    // ============================================================
    // Strategy: fill a 2D histogram (run, zDrift) via RDataFrame (MT-safe),
    // then after the event loop project each run bin, fit Gaussians, and
    // store results in TGraphErrors.

    // Define the run number as a column
    auto dfWithRun = dfZdrift.Define("RunNumber",
                                     [](ActRoot::MergerData& m) -> int { return m.fRun; },
                                     {"MergerData"});

    // 2D histo: x = run number, y = zDrift — bin width 1 in run so each bin = one run
    auto hZdriftVsRun = dfWithRun.Histo2D(
        {"hZdriftVsRun", "Z drift vs Run;Run number;Z drift [mm]", 300, 0, 300, 340, -40, 300},
        "RunNumber", "zDrift");

    // Force event loop
    hZdriftVsRun->GetEntries();

    // Containers for per-run fit results
    std::vector<double> runVec, meanStartVec, meanStartErrVec, meanEndVec, meanEndErrVec;
    std::vector<double> sigmaStartVec, sigmaStartErrVec, sigmaEndVec, sigmaEndErrVec;

    // Group runs in pairs of nGroup to increase statistics per fit point
    constexpr int nGroup = 1;
    int nBins = hZdriftVsRun->GetNbinsX();
    for(int ib = 1; ib <= nBins; ib += nGroup)
    {
        int ibEnd = std::min(ib + nGroup - 1, nBins);
        // Project the Y slice merging nGroup consecutive run bins
        TH1D* hSlice = hZdriftVsRun->ProjectionY(Form("hSlice_%d", ib), ib, ibEnd);
        if(hSlice->GetEntries() < 50) // need enough statistics
        {
            delete hSlice;
            continue;
        }

        // Run label = average of the bin-center range
        double runCenter = 0.5 * (hZdriftVsRun->GetXaxis()->GetBinCenter(ib) +
                                  hZdriftVsRun->GetXaxis()->GetBinCenter(ibEnd));

        // Fit pad-plane peak (low-Z extreme)
        TF1 fStart("fStart", "gaus", -20, 5);
        int statusStart = hSlice->Fit(&fStart, "RQ0");

        // Fit cathode peak (high-Z extreme)
        TF1 fEnd("fEnd", "gaus", 250, 270);
        int statusEnd = hSlice->Fit(&fEnd, "RQ0");

        // Only store if both fits converged
        if(statusStart == 0 && statusEnd == 0)
        {
            runVec.push_back(runCenter);

            meanStartVec.push_back(fStart.GetParameter(1));
            meanStartErrVec.push_back(fStart.GetParError(1));
            sigmaStartVec.push_back(fStart.GetParameter(2));
            sigmaStartErrVec.push_back(fStart.GetParError(2));

            meanEndVec.push_back(fEnd.GetParameter(1));
            meanEndErrVec.push_back(fEnd.GetParError(1));
            sigmaEndVec.push_back(fEnd.GetParameter(2));
            sigmaEndErrVec.push_back(fEnd.GetParError(2));
        }

        delete hSlice;
    }

    int nRuns = (int)runVec.size();
    std::vector<double> runErr(nRuns, 0.0); // no error on run number

    auto* grMeanStart = new TGraphErrors(nRuns, runVec.data(), meanStartVec.data(),
                                         runErr.data(), meanStartErrVec.data());
    grMeanStart->SetName("grMeanStart");
    grMeanStart->SetTitle("Pad-plane peak mean vs Run;Run number;#mu_{start} [mm]");
    grMeanStart->SetMarkerStyle(20);
    grMeanStart->SetMarkerColor(kRed);
    grMeanStart->SetLineColor(kRed);

    auto* grMeanEnd = new TGraphErrors(nRuns, runVec.data(), meanEndVec.data(),
                                       runErr.data(), meanEndErrVec.data());
    grMeanEnd->SetName("grMeanEnd");
    grMeanEnd->SetTitle("Cathode peak mean vs Run;Run number;#mu_{end} [mm]");
    grMeanEnd->SetMarkerStyle(20);
    grMeanEnd->SetMarkerColor(kBlue);
    grMeanEnd->SetLineColor(kBlue);

    auto* grSigmaStart = new TGraphErrors(nRuns, runVec.data(), sigmaStartVec.data(),
                                          runErr.data(), sigmaStartErrVec.data());
    grSigmaStart->SetName("grSigmaStart");
    grSigmaStart->SetTitle("Pad-plane peak #sigma vs Run;Run number;#sigma_{start} [mm]");
    grSigmaStart->SetMarkerStyle(21);
    grSigmaStart->SetMarkerColor(kRed);
    grSigmaStart->SetLineColor(kRed);

    auto* grSigmaEnd = new TGraphErrors(nRuns, runVec.data(), sigmaEndVec.data(),
                                        runErr.data(), sigmaEndErrVec.data());
    grSigmaEnd->SetName("grSigmaEnd");
    grSigmaEnd->SetTitle("Cathode peak #sigma vs Run;Run number;#sigma_{end} [mm]");
    grSigmaEnd->SetMarkerStyle(21);
    grSigmaEnd->SetMarkerColor(kBlue);
    grSigmaEnd->SetLineColor(kBlue);

    // Save zDrift per-run fit results to .dat
    {
        std::ofstream fout("./Outputs/zDrift_perRun.dat");
        fout << "# run  meanStart  sigmaStart  meanEnd  sigmaEnd\n";
        for(int i = 0; i < nRuns; ++i)
            fout << runVec[i] << "  " << meanStartVec[i] << "  " << sigmaStartVec[i]
                 << "  " << meanEndVec[i] << "  " << sigmaEndVec[i] << "\n";
        fout.close();
        std::cout << "Saved zDrift per-run fits to ./Outputs/zDrift_perRun.dat (" << nRuns << " entries)\n";
    }

    // Draw per-run fit evolution
    auto* cPerRun = new TCanvas("cPerRun", "zDrift extreme fits vs Run", 1200, 800);
    cPerRun->Divide(2, 2);

    cPerRun->cd(1);
    grMeanStart->Draw("APE");

    cPerRun->cd(2);
    grMeanEnd->Draw("APE");

    cPerRun->cd(3);
    grSigmaStart->Draw("APE");

    cPerRun->cd(4);
    grSigmaEnd->Draw("APE");

    // Also draw the 2D histo
    auto* cZdriftRun = new TCanvas("cZdriftRun", "Z drift vs Run (2D)", 800, 600);
    hZdriftVsRun->DrawClone("COLZ");

    // ============================================================
    // Per-run fit of Qave (same strategy as zDrift)
    // ============================================================
    // 2D histo: x = run number, y = Qave
    auto hQaveVsRun = dfWithRun.Histo2D(
        {"hQaveVsRun", "Qave vs Run;Run number;Qave [a.u.]", 300, 0, 300, 300, 0, 3000},
        "RunNumber", "MergerData.fQave");

    // Force event loop
    hQaveVsRun->GetEntries();

    // Containers for per-run Qave fit results
    std::vector<double> runQVec, meanQVec, meanQErrVec, sigmaQVec, sigmaQErrVec;

    for(int ib = 1; ib <= hQaveVsRun->GetNbinsX(); ib += nGroup)
    {
        int ibEnd = std::min(ib + nGroup - 1, hQaveVsRun->GetNbinsX());
        TH1D* hSliceQ = hQaveVsRun->ProjectionY(Form("hSliceQ_%d", ib), ib, ibEnd);
        if(hSliceQ->GetEntries() < 50)
        {
            delete hSliceQ;
            continue;
        }

        double runCenterQ = 0.5 * (hQaveVsRun->GetXaxis()->GetBinCenter(ib) +
                                   hQaveVsRun->GetXaxis()->GetBinCenter(ibEnd));

        // Fit Qave peak with Gaussian in the same range as the global fit
        TF1 fQ("fQ", "gaus", 0, 120);
        int statusQ = hSliceQ->Fit(&fQ, "RQ0");

        if(statusQ == 0)
        {
            runQVec.push_back(runCenterQ);
            meanQVec.push_back(fQ.GetParameter(1));
            meanQErrVec.push_back(fQ.GetParError(1));
            sigmaQVec.push_back(fQ.GetParameter(2));
            sigmaQErrVec.push_back(fQ.GetParError(2));
        }

        delete hSliceQ;
    }

    int nRunsQ = (int)runQVec.size();
    std::vector<double> runQErr(nRunsQ, 0.0);

    auto* grMeanQ = new TGraphErrors(nRunsQ, runQVec.data(), meanQVec.data(),
                                     runQErr.data(), meanQErrVec.data());
    grMeanQ->SetName("grMeanQ");
    grMeanQ->SetTitle("Qave peak mean vs Run;Run number;#mu_{Qave} [a.u.]");
    grMeanQ->SetMarkerStyle(20);
    grMeanQ->SetMarkerColor(kGreen + 2);
    grMeanQ->SetLineColor(kGreen + 2);

    auto* grSigmaQ = new TGraphErrors(nRunsQ, runQVec.data(), sigmaQVec.data(),
                                      runQErr.data(), sigmaQErrVec.data());
    grSigmaQ->SetName("grSigmaQ");
    grSigmaQ->SetTitle("Qave peak #sigma vs Run;Run number;#sigma_{Qave} [a.u.]");
    grSigmaQ->SetMarkerStyle(21);
    grSigmaQ->SetMarkerColor(kGreen + 2);
    grSigmaQ->SetLineColor(kGreen + 2);

    // Save Qave per-run fit results to .dat
    {
        std::ofstream fout("./Outputs/Qave_perRun.dat");
        fout << "# run  meanQave  sigmaQave\n";
        for(int i = 0; i < nRunsQ; ++i)
            fout << runQVec[i] << "  " << meanQVec[i] << "  " << sigmaQVec[i] << "\n";
        fout.close();
        std::cout << "Saved Qave per-run fits to ./Outputs/Qave_perRun.dat (" << nRunsQ << " entries)\n";
    }

    // Draw per-run Qave fit evolution
    auto* cPerRunQ = new TCanvas("cPerRunQ", "Qave fits vs Run", 1200, 400);
    cPerRunQ->Divide(2, 1);

    cPerRunQ->cd(1);
    grMeanQ->Draw("APE");

    cPerRunQ->cd(2);
    grSigmaQ->Draw("APE");

    // 2D Qave vs Run
    auto* cQaveRun = new TCanvas("cQaveRun", "Qave vs Run (2D)", 800, 600);
    hQaveVsRun->DrawClone("COLZ");

    // Save some specific data for debuging

    // cuts.ReadCut("lowQ", "./Cuts/events_lowQ_deposition.root");
    // dfFilter.Foreach(
    //     [&](ActRoot::MergerData& m)
    //     {
    //         if(cuts.IsInside("lowQ", m.fLight.fRawTL, m.fLight.fQtotal))
    //         {
    //             m.Stream(outFileCharge);
    //         }
    //     },
    //     {"MergerData"});
    // ActRoot::CutsManager<std::string> cuts;
    // cuts.ReadCut("lowQ", "./Cuts/events_lowQ_deposition.root");
    // std::ofstream outFileCharge("./Outputs/events_driftNearPad.dat");
    // dfZdrift.Foreach(
    //     [&](ActRoot::MergerData& m, float zDrift)
    //     {
    //         // if(cuts.IsInside("lowQ", m.fLight.fRawTL, m.fLight.fQtotal))
    //         // {
    //         if((zDrift > 0 && zDrift < 4)) // aprox pad plane and cathode position
    //         {
    //             m.Stream(outFileCharge);
    //         }
    //         // }
    //         return true;
    //     },
    //     {"MergerData", "zDrift"});
    // std::ofstream outFileCharge1("./Outputs/events_driftNearCathode.dat");
    // dfZdrift.Foreach(
    //     [&](ActRoot::MergerData& m, float zDrift)
    //     {
    //         // if(cuts.IsInside("lowQ", m.fLight.fRawTL, m.fLight.fQtotal))
    //         // {
    //         if((zDrift > 251 && zDrift < 254)) // aprox pad plane and cathode position
    //         {
    //             m.Stream(outFileCharge1);
    //         }
    //         // }
    //         return true;
    //     },
    //     {"MergerData", "zDrift"});
    // cuts.ReadCut("lowQ", "./Cuts/cut_lowQ_QaveVSZdrift.root");
    // std::ofstream outFileCharge2("./Outputs/events_lowQ.dat");
    // dfZdrift.Foreach(
    //     [&](ActRoot::MergerData& m, float zDrift)
    //     {
    //         if(cuts.IsInside("lowQ", m.fLight.fQave, zDrift))
    //         {
    //             m.Stream(outFileCharge2);
    //         }
    //         return true;
    //     },
    //     {"MergerData", "zDrift"});
}