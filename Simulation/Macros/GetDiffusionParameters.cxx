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
    // dataman.SetRuns(64, 67);
    // dataman.SetRuns(105, 122);
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
    auto dfIni = df.Filter(lambdaIsL1, {"MergerData", "ModularData"});

    auto lambdaIsSilLat {[](ActRoot::MergerData& mer, ActRoot::ModularData& mod)
                         { return (mod.Get("GATCONF") == 1 || mod.Get("GATCONF") == 2) && (mer.fLightIdx != -1); }};
    // auto dfIni = df.Filter(lambdaIsSilLat, {"MergerData", "ModularData"});

    // Filter events to avoid bad events for heavy or L1 analysis
    // Split filters into steps so we can count how many events survive each stage
    auto dfFilter1 = dfIni.Filter( // Check for heavier clusters than Li
        [](ActRoot::MergerData& m)
        {
            if(m.fHeavy.fQave > 1600.)
                return false;
            return true;
        },
        {"MergerData"});

    auto dfFilter2 = dfFilter1.Filter( // Check if heavy reaches end of ACTAR
        [](ActRoot::MergerData& m)
        {
            auto TL {m.fHeavy.fTL};
            auto theta {m.fThetaHeavy * TMath::DegToRad()};
            auto z_end {m.fRP.X() + TL * TMath::Cos(theta)};
            if(z_end < 240)
                return false;
            return true;
        },
        {"MergerData"});

    auto dfFilter3 = dfFilter2.Filter( // Mask events with large local charge near RP / beam cluster
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
                    if(v.GetPosition().Y() > rp_y - 3 && v.GetPosition().Y() < rp_y + 3) // aprox L1 exclusion zone
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
    auto dfFilterFinal = dfFilter3.Filter(
        [](ActRoot::MergerData& mer)
        {
            if(mer.fLight.fQave > 100.0 && mer.fLight.fQave < 1600.)
                return true;
            else
                return false;
        },
        {"MergerData"});

    // Print counts at each step to diagnose where events are lost
    auto nTotal = df.Count();
    auto nIni = dfIni.Count();
    auto nF1 = dfFilter1.Count();
    auto nF2 = dfFilter2.Count();
    auto nF3 = dfFilter3.Count();
    auto nFinal = dfFilterFinal.Count();

    std::cout << "Event counts:\n";
    std::cout << "  Total chain: " << *nTotal << "\n";
    std::cout << "  After GATCONF (1|2) & lightIdx!=-1: " << *nIni << "\n";
    std::cout << "  After fHeavy.fQave <= 1600: " << *nF1 << "\n";
    std::cout << "  After z_end >= 240: " << *nF2 << "\n";
    std::cout << "  After local high-charge mask near RP: " << *nF3 << "\n";
    std::cout << "  After fLight.fQave > 300: " << *nFinal << "\n";

    // Filter out events where the light-cluster fit is poor (keep chi2 <= 0.5)
    auto dfChi = dfFilterFinal.Filter(
        [](ActRoot::MergerData& mer, ActRoot::TPCData& tpc)
        {
            int idx = mer.fLightIdx;
            if(idx < 0)
                return false;
            const auto& cluster = tpc.fClusters.at(idx);
            double chi2 = cluster.GetLine().GetChi2();
            return (chi2 <= 0.4);
        },
        {"MergerData", "TPCData"});

    auto nAfterChi = dfChi.Count();
    std::cout << "  After chi2 <= 0.4 filter: " << *nAfterChi << "\n";

    // Get events that have high chi to visualize them
    // auto dfBadChi = dfFilterFinal.Filter(
    //     [](ActRoot::MergerData& mer, ActRoot::TPCData& tpc)
    //     {
    //         int idx = mer.fLightIdx;
    //         if(idx < 0)
    //             return false;
    //
    //         const auto& cluster = tpc.fClusters.at(idx);
    //         double chi2 = cluster.GetLine().GetChi2();
    //
    //         return (chi2 > 0.4);
    //     },
    //     {"MergerData", "TPCData"});

    // Dummy parameters. Need tuning
    int nBinsS = 15;
    double sMin = 0;
    double sMax = 300;
    int minVoxelsPerSlice = 8;
    double ds = (sMax - sMin) / nBinsS;

    auto dfSigmaLightPreFilter = dfChi.Define(
        "sigmaTransZ",
        [&](ActRoot::TPCData& tpc, ActRoot::MergerData& mer)
        {
            std::vector<std::pair<double, double>> out;

            int idx = mer.fLightIdx;
            if(idx < 0)
                return out;

            constexpr double zRPtoPad = 110.0; // mm (distance from RP to pad plane)

            auto cluster = tpc.fClusters.at(idx);
            auto line = cluster.GetLine();
            line.Scale(2, driftFactor); // scale line only (do NOT scale voxels)

            auto u = line.GetDirection().Unit();
            auto p0 = line.GetPoint(); // charge weighted center of the cluster

            for(int ib = 0; ib < nBinsS; ++ib)
            {
                double sLow = sMin + ib * ds;
                double sHigh = sLow + ds;

                double sumQ = 0.0;
                double sumR2 = 0.0;
                double sumZ_Q = 0.0;
                int nVox = 0;

                for(const auto& vox : cluster.GetVoxels())
                {
                    for(const auto& vext : vox.GetExtended())
                    {
                        // scaled position
                        const auto& pos = vext.GetPosition();
                        ROOT::Math::XYZPointF posScaled {pos.X() * 2.f, pos.Y() * 2.f, pos.Z() * driftFactor};

                        auto dr = posScaled - mer.fRP;
                        double s = dr.Dot(u);

                        if(s < sLow || s >= sHigh)
                            continue;

                        double q = vext.GetCharge();

                        auto d = posScaled - p0;
                        double par = d.Dot(u);
                        double r2 = d.Mag2() - par * par;

                        sumQ += q;
                        sumR2 += q * r2;
                        sumZ_Q += q * posScaled.Z();
                        nVox++;
                    }
                }

                // minimum amount of voxels and charge to consider the slice
                if(nVox < minVoxelsPerSlice || sumQ <= 0.)
                    continue;

                double sigma = std::sqrt(sumR2 / sumQ);

                // physical Z position of the slice (charge weighted)
                double zSlice = sumZ_Q / sumQ;

                // Z difference relative to the RP
                //  deltaZ > 0  -> track goes to larger Z (away from pad plane)
                //  deltaZ < 0  -> track goes to smaller Z (closer to pad plane)
                double deltaZ = zSlice - mer.fRP.Z();

                // real drift distance to the pad plane
                double zDrift = zRPtoPad + deltaZ;

                // optional physical protection
                // if(zDrift <= 0 || zDrift > zRPtoPad + 20)
                //     continue;

                out.emplace_back(zDrift, sigma);
            }

            return out;
        },
        {"TPCData", "MergerData"});

    // Filter particles that go very near pad plane or cathode (10 mm distance). They could have bad fit.
    auto dfSigmaLight = dfSigmaLightPreFilter.Filter(
        [](const std::vector<std::pair<double, double>>& v)
        {
            for(const auto& [z, sigma] : v)
            {
                if(z < 10 || z > 246)
                    return false;
            }
            return true;
        },
        {"sigmaTransZ"});

    // Print number of events that survive all filters and have sigmaTransZ info
    auto nSigmaLight = dfSigmaLight.Count();
    std::cout << "  After near-pad/cathode filter: " << *nSigmaLight << "\n";


    auto hSigmaZ = new TH2D("hSigmaZ", "#sigma_{trans} vs z (light);z from pad plane [mm];#sigma_{trans} [mm]", nBinsS,
                            sMin, sMax, 500, 0, 5);

    dfSigmaLight.Foreach(
        [&](const std::vector<std::pair<double, double>>& v)
        {
            for(const auto& [z, sigma] : v)
                hSigmaZ->Fill(z, sigma);
        },
        {"sigmaTransZ"});

    auto prof = hSigmaZ->ProfileX();

    // Do also the profile of sigma vs sqrt z and sigma² against z
    auto hSigmaSqrtZ = new TH2D(
        "hSigmaSqrtZ", "#sigma_{trans} vs sqrt(z) (light);sqrt(z) from pad plane [mm^{1/2}];#sigma_{trans} [mm]",
        std::sqrt(nBinsS), std::sqrt(sMin), std::sqrt(sMax), 500, 0, 5);
    auto hSigmaZ2 = new TH2D("hSigmaZ2", "#sigma^2_{trans} vs z (light);z from pad plane [mm];#sigma_{trans}^2 [mm^2]",
                             25, 0, 300, 500, 0, 25);
    dfSigmaLight.Foreach(
        [&](const std::vector<std::pair<double, double>>& v)
        {
            for(const auto& [z, sigma] : v)
            {
                hSigmaSqrtZ->Fill(std::sqrt(z), sigma);
                hSigmaZ2->Fill(z, sigma * sigma);
            }
        },
        {"sigmaTransZ"});

    // Get profiles and fit them lineary, then show results in canvas y = a + bx
    auto profSqrtZ = hSigmaSqrtZ->ProfileX();
    int sMinFitSqrtZ = 5;
    int sMaxFitSqrtZ = 13;
    profSqrtZ->Fit("pol1", "", "", sMinFitSqrtZ, sMaxFitSqrtZ);

    auto profZ2 = hSigmaZ2->ProfileX();
    // int sMinFitZ2 = 60;
    // int sMaxFitZ2 = 200;
    int sMinFitZ2 = 0;
    int sMaxFitZ2 = 250;
    profZ2->Fit("pol1", "", "", sMinFitZ2, sMaxFitZ2);

    // Debug z distances (I think they are too big)
    // auto dfDebug = dfSigmaLight
    //                    .Define("zBeginEnd_raw",
    //                            [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
    //                            {
    //                                int idx = m.fLightIdx;
    //                                if(idx < 0)
    //                                    return std::make_pair(0.f, 0.f);
    //
    //                                auto cluster = tpc.fClusters.at(idx);
    //                                auto line = cluster.GetLine();
    //                                auto u = line.GetDirection().Unit();
    //                                cluster.SortAlongDir(u);
    //                                float zMin = 1e6;
    //                                float zMax = -1e6;
    //                                for(auto& vox : cluster.GetVoxels())
    //                                {
    //                                    for(auto& vext : vox.GetExtended())
    //                                    {
    //                                        auto pos = vext.GetPosition();
    //                                        if(pos.Z() < zMin)
    //                                            zMin = pos.Z();
    //                                        if(pos.Z() > zMax)
    //                                            zMax = pos.Z();
    //                                    }
    //                                }
    //                                return std::make_pair(zMin * 4, zMax * 4); // Scale from btb to tb
    //                            },
    //                            {"MergerData", "TPCData"})
    //                    .Define("zMinMax_scaled",
    //                            [&](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
    //                            {
    //                                int idx = m.fLightIdx;
    //                                if(idx < 0)
    //                                    return std::make_pair(0.f, 0.f);
    //
    //                                auto cluster = tpc.fClusters.at(idx);
    //                                auto line = cluster.GetLine();
    //                                line.Scale(2, driftFactor); // scale line and not whole cluster to avoid moving
    //                                                            // voxels twice with GetExtended
    //                                auto u = line.GetDirection().Unit();
    //                                cluster.SortAlongDir(u);
    //                                float zMin = 1e6;
    //                                float zMax = -1e6;
    //                                for(auto& vox : cluster.GetVoxels())
    //                                {
    //                                    for(auto& vext : vox.GetExtended())
    //                                    {
    //                                        auto pos = vext.GetPosition();
    //                                        if(pos.Z() < zMin)
    //                                            zMin = pos.Z();
    //                                        if(pos.Z() > zMax)
    //                                            zMax = pos.Z();
    //                                    }
    //                                }
    //                                return std::make_pair(zMin, zMax);
    //                            },
    //                            {"MergerData", "TPCData"})
    //                    .Define("zDistance", [&](const std::pair<float, float>& zBeginEnd)
    //                            { return zBeginEnd.second - zBeginEnd.first; }, {"zMinMax_scaled"});
    //
    // // Diagnostic counts for zDistance to understand why many events may be missing from the
    // // hzDistance histogram (which is defined in the 0..300 mm range).
    // auto nDFDebug = dfDebug.Count();
    // auto nZD_inRange = dfDebug.Filter([](const float z) { return (z >= 0.0f && z <= 300.0f); },
    // {"zDistance"}).Count(); auto nZD_neg = dfDebug.Filter([](const float z) { return (z < 0.0f); },
    // {"zDistance"}).Count(); auto nZD_gt300 = dfDebug.Filter([](const float z) { return (z > 300.0f); },
    // {"zDistance"}).Count();
    //
    // std::cout << "dfDebug counts:\n";
    // std::cout << "  Total dfDebug rows: " << *nDFDebug << "\n";
    // std::cout << "  zDistance in [0,300]: " << *nZD_inRange << "\n";
    // std::cout << "  zDistance < 0: " << *nZD_neg << "\n";
    // std::cout << "  zDistance > 300: " << *nZD_gt300 << "\n";

    auto hzBeginRaw = new TH1D("hzBeginRaw", "z begin (raw);z [mm]", 100, 0, 512);
    auto hzEndRaw = new TH1D("hzEndRaw", "z end (raw);z [mm]", 100, 0, 512);
    auto hzBeginScaled = new TH1D("hzBeginScaled", "z begin (scaled);z [mm]", 100, 0, 512);
    auto hzEndScaled = new TH1D("hzEndScaled", "z end (scaled);z [mm]", 100, 0, 512);
    auto hzDistance = new TH1D("hzDistance", "z distance (scaled);z [mm]", 100, 0, 300);

    // dfDebug.Foreach(
    //     [&](const std::pair<float, float>& zBeginEnd, const std::pair<float, float>& zMinMaxScaled,
    //         const float zDistance)
    //     {
    //         hzBeginScaled->Fill(zMinMaxScaled.first);
    //         hzEndScaled->Fill(zMinMaxScaled.second);
    //
    //        hzBeginRaw->Fill(zBeginEnd.first);
    //        hzEndRaw->Fill(zBeginEnd.second);
    //
    //        hzDistance->Fill(zDistance);
    //    },
    //    {"zBeginEnd_raw", "zMinMax_scaled", "zDistance"});

    // Debug strange behaviour at distZ < than 110 mm
    // int causalBreaks {0};
    // dfSigmaLight.Foreach(
    //    [&](const std::vector<std::pair<double, double>>& slices, ActRoot::MergerData& mer)
    //    {
    //        constexpr double zRPtoPad = 110.0; // same as in Define
    //
    //        for(const auto& [zDrift, sigma] : slices)
    //        {
    //            // deltaZ = zSlice - RP.Z() = zDrift - zRPtoPad
    //            double deltaZ = zDrift - zRPtoPad;
    //
    //            // check if deltaZ positive but slice too close to pad plane
    //            if(deltaZ > 0 && zDrift < 110.0)
    //                causalBreaks++;
    //        }
    //    },
    //    {"sigmaTransZ", "MergerData"});
    //
    // std::cout << "Number of slices breaking causality (deltaZ > 0 but zDrift < 110 mm): " << causalBreaks <<
    // std::endl;

    // Get histos for zdrift distribution and phi angle distribution to check posible assymetry
    auto hZdrift = new TH1D("hZdrift", "Drift distance to pad plane;z_{drift} [mm];Entries", 100, 0, 300);
    auto hPhiLight = new TH1D("hPhiLight", "Phi angle of light particle;#phi_{light} [deg];Entries", 50, -180, 180);
    auto hThetaLight =
        new TH1D("hThetaLight", "Theta angle of light particle;#theta_{light} [deg];Entries", 50, 0, 180);
    auto hDirZ = new TH1D("hDirZ", "Direction Z component;u_{z};Counts", 200, -1, 1);
    dfSigmaLight.Foreach(
        [&](const std::vector<std::pair<double, double>>& slices, ActRoot::MergerData& mer, ActRoot::TPCData& tpc)
        {
            constexpr double zRPtoPad = 110.0; // same as in Define

            // Get only the last slice (final position of particle)
            if(slices.empty())
                return;

            const auto& [zDrift, sigma] = slices.back();
            hZdrift->Fill(zDrift);

            int idx = mer.fLightIdx;
            auto line = tpc.fClusters[idx].GetLine();
            auto u = line.GetDirection().Unit();
            hDirZ->Fill(u.Z());

            hPhiLight->Fill(mer.fPhiLight);

            hThetaLight->Fill(mer.fThetaLight);
        },
        {"sigmaTransZ", "MergerData", "TPCData"});

    // Get slices number for particles in silcon (a lot of them empty I think)
    auto hNSlices = new TH1D("hNSlices", "Number of valid slices per event;N slices;Events", 20, 0, 20);

    dfSigmaLight.Foreach([&](const std::vector<std::pair<double, double>>& v) { hNSlices->Fill(v.size()); },
                         {"sigmaTransZ"});


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

    gStyle->SetOptFit(1111);
    auto* c1 = new TCanvas("c1", "Transverse sigma of light particle", 800, 600);
    c1->Divide(2, 1);
    c1->cd(1);
    hSigmaZ->DrawClone();
    c1->cd(2);
    prof->DrawClone();

    auto* c2 = new TCanvas("c2", "Transverse sigma of light particle FIT", 800, 600);
    c2->Divide(2, 2);
    c2->cd(1);
    hSigmaSqrtZ->DrawClone();
    c2->cd(2);
    profSqrtZ->DrawClone();
    c2->cd(3);
    hSigmaZ2->DrawClone();
    c2->cd(4);
    profZ2->DrawClone();

    auto* c3 = new TCanvas("c3", "Z distances", 800, 600);
    c3->Divide(2, 2);
    c3->cd(1);
    hzBeginRaw->DrawClone();
    c3->cd(2);
    hzEndRaw->DrawClone();
    c3->cd(3);
    hzBeginScaled->DrawClone();
    c3->cd(4);
    hzEndScaled->DrawClone();

    auto* c4 = new TCanvas("c4", "Z distance between begin and end", 800, 600);
    hzDistance->DrawClone();

    auto* c5 = new TCanvas("c5", "Z drift and phi distribution", 800, 600);
    c5->Divide(2, 2);
    c5->cd(1);
    hZdrift->DrawClone();
    c5->cd(2);
    hPhiLight->DrawClone();
    c5->cd(3);
    hThetaLight->DrawClone();
    c5->cd(4);
    hDirZ->DrawClone();

    auto* c6 = new TCanvas("c6", "Number of slices", 800, 600);
    hNSlices->DrawClone();

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
    // Save events with too much distance in z
    // std::ofstream outFile1 {"./Outputs/Events_HighZ.dat"};
    // dfSigmaLight.Foreach(
    //     [&outFile1](const ActRoot::MergerData& mer, const std::vector<std::pair<double, double>>& v)
    //     {
    //     for(const auto& [z, sigma] : v)
    //     {
    //         if(z > 300)
    //         {
    //             mer.Stream(outFile1);
    //             break; // Only need to save the event once, even if it has multiple slices with s < 0
    //         }
    //     }
    //     },
    //     {"MergerData", "sigmaTransZ"});
    // outFile1.close();
    std::ofstream outFile2 {"./Outputs/Events_DriftDistLower60SigmaHigher1-4_NoNearBorders_L1.dat"};
    dfSigmaLight.Foreach(
        [&outFile2](const ActRoot::MergerData& mer, const std::vector<std::pair<double, double>>& v)
        {
            for(const auto& [zDistance, sigma] : v)
            {
                if((zDistance < 70 && sigma > 1.4) || (zDistance < 130 && sigma > 1.3))
                {
                    mer.Stream(outFile2);
                    break; // Only need to save the event once, even if it has multiple slices with s < 0
                }
            }
        },
        {"MergerData", "sigmaTransZ"});
    outFile2.close();
    std::ofstream outFile3 {"./Outputs/Events_DriftDistGreater190SigmaGreater2_NoNearBorders_L1.dat"};
    dfSigmaLight.Foreach(
        [&outFile3](const ActRoot::MergerData& mer, const std::vector<std::pair<double, double>>& v)
        {
            for(const auto& [zDistance, sigma] : v)
            {
                if(((zDistance > 200 && zDistance < 230) && sigma > 1.8) ||
                   ((zDistance > 140 && zDistance < 160) && sigma > 1.6))
                {
                    mer.Stream(outFile3);
                    break; // Only need to save the event once, even if it has multiple slices with s < 0
                }
            }
        },
        {"MergerData", "sigmaTransZ"});
    // outFile3.close();
    // std::ofstream outFile4 {"./Outputs/Events_ChiGreater0-4_L1beforeThresholdChange.dat"};
    // dfBadChi.Foreach([&outFile4](const ActRoot::MergerData& mer, const ActRoot::TPCData& tpc) { mer.Stream(outFile4);
    // },
    //                  {"MergerData", "TPCData"});
}