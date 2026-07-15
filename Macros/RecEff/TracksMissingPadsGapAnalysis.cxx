// TracksMissingPadsGapAnalysis.C
//
// Studies gaps along the s-coordinate obtained by projecting the light
// cluster's voxels onto its fitted line. Also defines the good/bad event
// classification (via CutsManager on maxGap vs meanGap) reused for the PID,
// theta and phi plots.
//
// This is one half of what used to be a single macro; the XZ/YZ pad
// occupancy study now lives in TracksMissingPadsXZYZProjection.C and shares
// the struct/function definitions in GapAnalysisHelpers.h.

#include "MissingPads.h"

#include "ActContinuity.h"
#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActLine.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActTPCData.h"
#include "ActTPCParameters.h"
#include "ActVoxel.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TROOT.h"

#include "Math/Point3Dfwd.h"
#include "Math/Vector3D.h"

#include <fstream>
#include <iostream>
#include <string>

void TracksMissingPadsGapAnalysis()
{
    ROOT::EnableImplicitMT(); // Multithreading is worthwhile since we process the full dataset.

    // ActRoot::DataManager dataman {"../../configs/data_7Li.conf", ActRoot::ModeType::EMerge};
    // dataman.SetRuns(69, 69);
    // auto chain {dataman.GetChain()};
    // auto friend1 {dataman.GetChain(ActRoot::ModeType::EFilter)};
    // chain->AddFriend(friend1.get());
    // auto friend2 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    // chain->AddFriend(friend2.get());

    ActRoot::InputParser parser {};
    parser.ReadFile("../../configs/detector.conf");
    auto driftBlock = parser.GetBlock("Merger");
    auto driftFactor = driftBlock->GetDouble("DriftFactor");

    // ROOT::RDataFrame df {*chain};

    ROOT::RDataFrame df {"PreProcessed_Tree", "../../PostAnalysis/Outputs/tree_preprocess_F_11Li.root"};

    auto dfFiltered = df.Filter([](ActRoot::ModularData& m, ActRoot::MergerData& merger)
                                { return m.Get("GATCONF") == 8; }, {"ModularData", "MergerData"}); // only L1

    // Previously done with && merger.fRun == 69

    // No event selection: process the entire chain.
    // If needed, any physics selection (e.g. only L1 events) can be added here.
    auto defSummary =
        dfFiltered
            .Define("projPoints",
                    [&driftFactor](ActRoot::MergerData& m, ActRoot::TPCData tpc)
                    {
                        std::vector<ProjectedPoint> empty;
                        auto lightIdx = m.fLightIdx;
                        if(lightIdx < 0 || lightIdx >= (int)tpc.fClusters.size())
                            return empty;
                        auto lightCl = tpc.fClusters[lightIdx];
                        auto rp = m.fRP;
                        return GetProjectedPositionsMM(lightCl, rp, driftFactor);
                    },
                    {"MergerData", "TPCData"})
            .Define("gapSummary", [](const std::vector<ProjectedPoint>& points) { return ComputeGapSummary(points); },
                    {"projPoints"})
            .Define("maxGap", [](const GapSummary& g) { return g.maxGap; }, {"gapSummary"})
            .Define("maxGapNorm", [](const GapSummary& g) { return g.maxGapNorm; }, {"gapSummary"})
            .Define("medianGap", [](const GapSummary& g) { return g.medianGap; }, {"gapSummary"})
            .Define("meanGap", [](const GapSummary& g) { return g.meanGap; }, {"gapSummary"})
            .Define("nVoxelsG", [](const GapSummary& g) { return g.nVoxels; }, {"gapSummary"})
            .Define("sRange", [](const GapSummary& g) { return g.sRange; }, {"gapSummary"})
            .Define("stdGap", [](const GapSummary& g) { return g.stdGap; }, {"gapSummary"})
            .Define("mad", [](const GapSummary& g) { return g.mad; }, {"gapSummary"})
            .Define("maxGapZ", [](const GapSummary& g) { return g.maxGapZ; }, {"gapSummary"})
            .Define("skewness", [](const GapSummary& g) { return g.skewness; }, {"gapSummary"})
            .Filter([](const GapSummary& g) { return g.valid; },
                    {"gapSummary"}); // Discard events with fewer than 3 voxels.

    // Get the data in the two branches of phi
    auto defSummaryPhiPositive =
        defSummary.Filter([](ActRoot::MergerData& m) { return m.fPhiLight > 0; }, {"MergerData"});
    auto defSummaryPhiNegative =
        defSummary.Filter([](ActRoot::MergerData& m) { return m.fPhiLight < 0; }, {"MergerData"});
    auto defSummaryPhiSmaller90 =
        defSummary.Filter([](ActRoot::MergerData& m) { return std::abs(m.fPhiLight) < 90; }, {"MergerData"});
    auto defSummaryPhiLarger90 =
        defSummary.Filter([](ActRoot::MergerData& m) { return std::abs(m.fPhiLight) > 90; }, {"MergerData"});

    // --- Distribution of the absolute maximum gap (mm) ---
    auto hMaxGap = defSummary.Histo1D({"hMaxGap", "Max gap per event;max gap [mm];Events", 200, 0, 10}, "maxGap");
    auto hMaxGapPhiPositive = defSummaryPhiPositive.Histo1D(
        {"hMaxGapPhiPositive", "Max gap per event (phi>0);max gap [mm];Events", 200, 0, 10}, "maxGap");
    auto hMaxGapPhiNegative = defSummaryPhiNegative.Histo1D(
        {"hMaxGapPhiNegative", "Max gap per event (phi<0);max gap [mm];Events", 200, 0, 10}, "maxGap");

    // --- Distribution of the mean gap (mm) ---
    auto hMeanGap = defSummary.Histo1D({"hMeanGap", "Mean gap per event;mean gap [mm];Events", 200, 0, 3.5}, "meanGap");
    auto hMeanGapPhiPositive = defSummaryPhiPositive.Histo1D(
        {"hMeanGapPhiPositive", "Mean gap per event (phi>0);mean gap [mm];Events", 200, 0, 3.5}, "meanGap");
    auto hMeanGapPhiNegative = defSummaryPhiNegative.Histo1D(
        {"hMeanGapPhiNegative", "Mean gap per event (phi<0);mean gap [mm];Events", 200, 0, 3.5}, "meanGap");
    auto hMeanGapVSMaxGap = defSummary.Histo2D(
        {"hMeanGapVSMaxGap", "Mean gap vs max gap;max gap [mm];mean gap [mm]", 200, 0, 10, 200, 0, 3.5}, "maxGap",
        "meanGap");
    auto hMeanGapVSMaxGapPhiPositive = defSummaryPhiPositive.Histo2D(
        {"hMeanGapVSMaxGapPhiPositive", "Mean gap vs max gap (phi>0);max gap [mm];mean gap [mm]", 200, 0, 10, 200, 0,
         3.5},
        "maxGap", "meanGap");
    auto hMeanGapVSMaxGapPhiNegative = defSummaryPhiNegative.Histo2D(
        {"hMeanGapVSMaxGapPhiNegative", "Mean gap vs max gap (phi<0);max gap [mm];mean gap [mm]", 200, 0, 10, 200, 0,
         3.5},
        "maxGap", "meanGap");
    auto hMeanGapVSMaxGapPhiSmaller90 = defSummaryPhiSmaller90.Histo2D(
        {"hMeanGapVSMaxGapPhiSmaller90", "Mean gap vs max gap (|phi|<90);max gap [mm];mean gap [mm]", 200, 0, 10, 200,
         0, 3.5},
        "maxGap", "meanGap");
    auto hMeanGapVSMaxGapPhiLarger90 = defSummaryPhiLarger90.Histo2D(
        {"hMeanGapVSMaxGapPhiLarger90", "Mean gap vs max gap (|phi|>90);max gap [mm];mean gap [mm]", 200, 0, 10, 200, 0,
         3.5},
        "maxGap", "meanGap");

    // --- Distribution of the standard deviation of gaps (mm) ---
    auto hStdGap =
        defSummary.Histo1D({"hStdGap", "Standard deviation of gaps;std gap [mm];Events", 200, 0, 3.5}, "stdGap");
    auto hStdGapMaxGap = defSummary.Histo2D(
        {"hStdGapMaxGap", "Standard deviation of gaps vs max gap;max gap [mm];std gap [mm]", 200, 0, 10, 200, 0, 3.5},
        "maxGap", "stdGap");
    auto hStdGapMeanGap =
        defSummary.Histo2D({"hStdGapMeanGap", "Standard deviation of gaps vs mean gap;mean gap [mm];std gap [mm]", 200,
                            0, 3.5, 200, 0, 3.5},
                           "meanGap", "stdGap");
    auto hStdGapMedianGap =
        defSummary.Histo2D({"hStdGapMedianGap", "Standard deviation of gaps vs median gap;median gap [mm];std gap [mm]",
                            200, 0, 3.5, 200, 0, 3.5},
                           "medianGap", "stdGap");

    // --- Distribution of the robust z-score of maxGap (MAD-based) ---
    auto hMaxGapZ = defSummary.Histo1D({"hMaxGapZ", "Robust z-score of max gap;maxGapZ;Events", 200, 0, 50}, "maxGapZ");
    auto hMaxGapZvsMaxGap = defSummary.Histo2D(
        {"hMaxGapZvsMaxGap", "maxGapZ vs max gap;max gap [mm];maxGapZ", 200, 0, 10, 800, 0, 200}, "maxGap", "maxGapZ");
    auto hMaxGapZvsMeanGap =
        defSummary.Histo2D({"hMaxGapZvsMeanGap", "maxGapZ vs mean gap;mean gap [mm];maxGapZ", 200, 0, 3.5, 200, 0, 50},
                           "meanGap", "maxGapZ");
    auto hMaxGapZvsMedianGap = defSummary.Histo2D(
        {"hMaxGapZvsMedianGap", "maxGapZ vs median gap;median gap [mm];maxGapZ", 200, 0, 3.5, 200, 0, 50}, "medianGap",
        "maxGapZ");

    // --- Distribution of the skewness of gaps ---
    auto hSkewness = defSummary.Histo1D({"hSkewness", "Skewness of gaps;skewness;Events", 200, -5, 30}, "skewness");
    auto hSkewnessVsMaxGap =
        defSummary.Histo2D({"hSkewnessVsMaxGap", "Skewness vs max gap;max gap [mm];skewness", 200, 0, 10, 200, -5, 30},
                           "maxGap", "skewness");
    auto hSkewnessVsMeanGap = defSummary.Histo2D(
        {"hSkewnessVsMeanGap", "Skewness vs mean gap;mean gap [mm];skewness", 200, 0, 3.5, 200, -5, 30}, "meanGap",
        "skewness");
    auto hSkewnessVsMedianGap = defSummary.Histo2D(
        {"hSkewnessVsMedianGap", "Skewness vs median gap;median gap [mm];skewness", 200, 0, 3.5, 200, -5, 30},
        "medianGap", "skewness");

    // --- Distribution of the normalized maximum gap (dimensionless) ---
    // This is likely the most useful observable for defining a threshold.
    // A maxGapNorm around 1-2 corresponds to the expected spacing, whereas
    // values much larger than 1 indicate a suspicious discontinuity,
    // regardless of track length, orientation, or voxel density.
    auto hMaxGapNorm =
        defSummary.Histo1D({"hMaxGapNorm", "Max gap / median gap;maxGap / medianGap;Events", 200, 0, 50}, "maxGapNorm");

    // --- Correlation: does the maximum gap depend on the number of voxels? ---
    // This helps determine whether a global threshold is sufficient or whether
    // it should depend on the track size.
    auto hGapVsNVoxels = defSummary.Histo2D(
        {"hGapVsNVoxels", "Max gap vs nVoxels;nVoxels;max gap [mm]", 100, 0, 500, 100, 0, 50}, "nVoxelsG", "maxGap");

    auto hGapNormVsNVoxels = defSummary.Histo2D(
        {"hGapNormVsNVoxels", "Max gap norm vs nVoxels;nVoxels;maxGap/medianGap", 100, 0, 500, 100, 0, 50}, "nVoxelsG",
        "maxGapNorm");

    // --- Correlation with the total projected track length (sRange) ---
    // Longer tracks naturally have more opportunities to contain a large gap
    // purely by chance.
    auto hGapVsSRange = defSummary.Histo2D(
        {"hGapVsSRange", "Max gap vs track length;s range [mm];max gap [mm]", 100, 0, 300, 100, 0, 50}, "sRange",
        "maxGap");

    auto* c1 = new TCanvas("c1", "Max gap distributions", 1200, 800);
    c1->Divide(2, 2);
    c1->cd(1);
    hMaxGap->DrawClone();
    c1->cd(2);
    hMaxGapNorm->DrawClone();
    c1->cd(3);
    hGapVsNVoxels->DrawClone("colz");
    c1->cd(4);
    hGapNormVsNVoxels->DrawClone("colz");

    auto* c2 = new TCanvas("c2", "Max gap vs track length", 800, 600);
    hGapVsSRange->DrawClone("colz");

    auto* c3 = new TCanvas("c3", "Max gap distributions by phi", 1200, 400);
    c3->Divide(2, 1);
    c3->cd(1);
    hMaxGapPhiPositive->DrawClone();
    c3->cd(2);
    hMaxGapPhiNegative->DrawClone();

    auto* c33 = new TCanvas("c33", "Mean gap distributions", 1200, 400);
    c33->Divide(3, 2);
    c33->cd(1);
    hMeanGap->DrawClone();
    c33->cd(2);
    hMeanGapVSMaxGap->DrawClone("colz");
    c33->cd(3);
    hMeanGapVSMaxGapPhiPositive->DrawClone("colz");
    c33->cd(4);
    hMeanGapVSMaxGapPhiNegative->DrawClone("colz");
    c33->cd(5);
    hMeanGapVSMaxGapPhiSmaller90->DrawClone("colz");
    c33->cd(6);
    hMeanGapVSMaxGapPhiLarger90->DrawClone("colz");

    auto* c34 = new TCanvas("c34", "Standard deviation of gaps", 1200, 400);
    c34->Divide(2, 2);
    c34->cd(1);
    hStdGap->DrawClone();
    c34->cd(2);
    hStdGapMaxGap->DrawClone("colz");
    c34->cd(3);
    hStdGapMeanGap->DrawClone("colz");
    c34->cd(4);
    hStdGapMedianGap->DrawClone("colz");

    auto* c35 = new TCanvas("c35", "Robust z-score of max gap (MAD-based)", 1200, 400);
    c35->Divide(2, 2);
    c35->cd(1);
    hMaxGapZ->DrawClone();
    c35->cd(2);
    hMaxGapZvsMaxGap->DrawClone("colz");
    c35->cd(3);
    hMaxGapZvsMeanGap->DrawClone("colz");
    c35->cd(4);
    hMaxGapZvsMedianGap->DrawClone("colz");

    auto* c36 = new TCanvas("c36", "Skewness of gaps", 1200, 400);
    c36->Divide(2, 2);
    c36->cd(1);
    hSkewness->DrawClone();
    c36->cd(2);
    hSkewnessVsMaxGap->DrawClone("colz");
    c36->cd(3);
    hSkewnessVsMeanGap->DrawClone("colz");
    c36->cd(4);
    hSkewnessVsMedianGap->DrawClone("colz");

    std::cout << "Total events processed: " << *defSummary.Count() << std::endl;

    // Create cuts and save events inside
    ActRoot::CutsManager<std::string> cuts {};
    cuts.ReadCut("cut", "./Cuts/cut_LimitTest_maxGap_meanGap.root");
    c33->cd(2);
    cuts.DrawCut("cut");
    std::ofstream outFile("./Outputs/eventsMaxGapMeanGap_LimitTestEvents.dat");
    // defSummary.Foreach(
    //     [&outFile, &cuts](ActRoot::MergerData& m, double meanGap, double maxGap)
    //     {
    //         if(cuts.IsInside("cut", maxGap, meanGap))
    //         {
    //             m.Stream(outFile);
    //         }
    //     },
    //     {"MergerData", "meanGap", "maxGap"});

    // Save some events if needed
    // std::ofstream outFile("./Outputs/eventsMaxGap_More8.dat");
    // defSummary.Filter([](double maxGap) { return maxGap > 8.0; }, {"maxGap"})
    //     .Foreach([&outFile](ActRoot::MergerData& m, ActRoot::TPCData& tpc) { m.Stream(outFile); },
    //              {"MergerData", "TPCData"});

    // Do PID for the good events in cut
    cuts.ReadCut("goodEvents", "./Cuts/cut_goodEvents_maxGap_meanGap.root");
    cuts.ReadCut("goodEventsWider", "./Cuts/cut_goodEventsWider_maxGap_meanGap.root");
    cuts.ReadCut("notGoodEvents", "./Cuts/cut_NotGoodEvents_maxGap_meanGap.root");
    cuts.ReadCut("notAtAllGoodEvents", "./Cuts/cut_NotAtAllGoodEvents_maxGap_meanGap.root");
    cuts.DrawCut("goodEventsWider");
    cuts.DrawCut("goodEvents");
    cuts.DrawCut("notGoodEvents");
    cuts.DrawCut("notAtAllGoodEvents");
    auto dfGoodEvents =
        defSummary.Filter([&cuts](double meanGap, double maxGap)
                          { return cuts.IsInside("goodEvents", maxGap, meanGap); }, {"meanGap", "maxGap"});
    auto dfGoodEventsWider =
        defSummary.Filter([&cuts](double meanGap, double maxGap)
                          { return cuts.IsInside("goodEventsWider", maxGap, meanGap); }, {"meanGap", "maxGap"});
    auto dfNotGoodEvents =
        defSummary.Filter([&cuts](double meanGap, double maxGap)
                          { return cuts.IsInside("notGoodEvents", maxGap, meanGap); }, {"meanGap", "maxGap"});
    auto dfNotAtAllGoodEvents =
        defSummary.Filter([&cuts](double meanGap, double maxGap)
                          { return cuts.IsInside("notAtAllGoodEvents", maxGap, meanGap); }, {"meanGap", "maxGap"});

    // PID in Qtot vs TLraw
    auto hPIDgoodEvents = dfGoodEvents.Histo2D({"hPID", "PID for good events;TL;Qtot", 200, 0, 120, 2000, 0, 3e5},
                                               "fLight.fRawTL", "fLight.fQtotal");
    auto hPIDgoodEventsWider =
        dfGoodEventsWider.Histo2D({"hPIDwider", "PID for good events (wider cut);TL;Qtot", 200, 0, 120, 2000, 0, 3e5},
                                  "fLight.fRawTL", "fLight.fQtotal");
    auto hPIDnotGoodEvents =
        dfNotGoodEvents.Histo2D({"hPIDnotGood", "PID for not good events;TL;Qtot", 200, 0, 120, 2000, 0, 3e5},
                                "fLight.fRawTL", "fLight.fQtotal");
    auto hPIDnotAtAllGoodEvents = dfNotAtAllGoodEvents.Histo2D(
        {"hPIDnotAtAllGood", "PID for not at all good events;TL;Qtot", 200, 0, 120, 2000, 0, 3e5}, "fLight.fRawTL",
        "fLight.fQtotal");
    auto hPIDgoodEventsPhiPositive =
        dfGoodEvents.Filter([](ActRoot::MergerData& m) { return m.fPhiLight > 0; }, {"MergerData"})
            .Histo2D({"hPIDgoodEventsPhiPositive", "PID for good events (phi>0);TL;Qtot", 200, 0, 120, 2000, 0, 3e5},
                     "fLight.fRawTL", "fLight.fQtotal");
    auto hPIDgoodEventsPhiNegative =
        dfGoodEvents.Filter([](ActRoot::MergerData& m) { return m.fPhiLight < 0; }, {"MergerData"})
            .Histo2D({"hPIDgoodEventsPhiNegative", "PID for good events (phi<0);TL;Qtot", 200, 0, 120, 2000, 0, 3e5},
                     "fLight.fRawTL", "fLight.fQtotal");
    auto hPIDall = defSummary.Histo2D({"hPIDall", "PID for all events;TL;Qtot", 200, 0, 120, 2000, 0, 3e5},
                                      "fLight.fRawTL", "fLight.fQtotal");
    auto hPIDallPhiPositive =
        defSummary.Filter([](ActRoot::MergerData& m) { return m.fPhiLight > 0; }, {"MergerData"})
            .Histo2D({"hPIDallPhiPositive", "PID for all events (phi>0);TL;Qtot", 200, 0, 120, 2000, 0, 3e5},
                     "fLight.fRawTL", "fLight.fQtotal");
    auto hPIDallPhiNegative =
        defSummary.Filter([](ActRoot::MergerData& m) { return m.fPhiLight < 0; }, {"MergerData"})
            .Histo2D({"hPIDallPhiNegative", "PID for all events (phi<0);TL;Qtot", 200, 0, 120, 2000, 0, 3e5},
                     "fLight.fRawTL", "fLight.fQtotal");

    // PIDs in Qave vs TL
    auto hPIDgoodEventsQave = dfGoodEvents.Histo2D(
        {"hPIDQave", "PID for good events;TL;Qave", 500, 0, 300, 2000, 0, 3e3}, "fLight.fTL", "fLight.fQave");
    auto hPIDgoodEventsWiderQave = dfGoodEventsWider.Histo2D(
        {"hPIDwiderQave", "PID for good events (wider cut);TL;Qave", 500, 0, 300, 2000, 0, 3e3}, "fLight.fTL",
        "fLight.fQave");
    auto hPIDnotGoodEventsQave =
        dfNotGoodEvents.Histo2D({"hPIDnotGoodQave", "PID for not good events;TL;Qave", 500, 0, 300, 2000, 0, 3e3},
                                "fLight.fTL", "fLight.fQave");
    auto hPIDnotAtAllGoodEventsQave = dfNotAtAllGoodEvents.Histo2D(
        {"hPIDnotAtAllGoodQave", "PID for not at all good events;TL;Qave", 500, 0, 300, 2000, 0, 3e3}, "fLight.fTL",
        "fLight.fQave");

    // Theta and phi distribution for good, bad and all events
    auto hThetaGoodEvents = dfGoodEvents.Histo1D(
        {"hThetaGoodEvents", "Theta distribution for good events;theta [deg];Events", 180, 0, 180}, "fThetaLight");
    auto hThetaAll = defSummary.Histo1D(
        {"hThetaAll", "Theta distribution for all events;theta [deg];Events", 180, 0, 180}, "fThetaLight");
    auto hThetaBad =
        defSummary
            .Filter([&cuts](ActRoot::MergerData& m, double maxGap, double meanGap)
                    { return !cuts.IsInside("goodEvents", maxGap, meanGap); }, {"MergerData", "maxGap", "meanGap"})
            .Histo1D({"hThetaBad", "Theta distribution for bad events;theta [deg];Events", 180, 0, 180}, "fThetaLight");
    auto hPhiGoodEvents = dfGoodEvents.Histo1D(
        {"hPhiGoodEvents", "Phi distribution for good events;phi [deg];Events", 360, -180, 180}, "fPhiLight");
    auto hPhiAll = defSummary.Histo1D({"hPhiAll", "Phi distribution for all events;phi [deg];Events", 360, -180, 180},
                                      "fPhiLight");
    auto hPhiBad =
        defSummary
            .Filter([&cuts](ActRoot::MergerData& m, double maxGap, double meanGap)
                    { return !cuts.IsInside("goodEvents", maxGap, meanGap); }, {"MergerData", "maxGap", "meanGap"})
            .Histo1D({"hPhiBad", "Phi distribution for bad events;phi [deg];Events", 360, -180, 180}, "fPhiLight");

    auto* cPIDgoodAndAll = new TCanvas("cPID", "PID for good events", 800, 600);
    cPIDgoodAndAll->Divide(2, 1);
    cPIDgoodAndAll->cd(1);
    hPIDgoodEvents->DrawClone("colz");
    cPIDgoodAndAll->cd(2);
    hPIDall->DrawClone("colz");

    auto* cPIDphi = new TCanvas("cPIDphi", "PID for all events", 800, 600);
    cPIDphi->Divide(2, 2);
    cPIDphi->cd(1);
    hPIDgoodEventsPhiPositive->DrawClone("colz");
    cPIDphi->cd(2);
    hPIDgoodEventsPhiNegative->DrawClone("colz");
    cPIDphi->cd(3);
    hPIDallPhiPositive->DrawClone("colz");
    cPIDphi->cd(4);
    hPIDallPhiNegative->DrawClone("colz");

    auto* cPIDWiderComparison = new TCanvas("cPIDwider", "PID for good events (wider cut)", 800, 600);
    cPIDWiderComparison->Divide(2, 2);
    cPIDWiderComparison->cd(1);
    hPIDgoodEventsWider->DrawClone("colz");
    cPIDWiderComparison->cd(2);
    hPIDgoodEvents->DrawClone();
    cPIDWiderComparison->cd(3);
    hPIDnotGoodEvents->DrawClone("colz");
    cPIDWiderComparison->cd(4);
    hPIDnotAtAllGoodEvents->DrawClone("colz");

    auto* cPIDsPhysicalUnits = new TCanvas("cPIDsPhysicalUnits", "PID for good events (physical units)", 800, 600);
    cPIDsPhysicalUnits->Divide(2, 2);
    cPIDsPhysicalUnits->cd(1);
    hPIDgoodEventsQave->DrawClone("colz");
    cPIDsPhysicalUnits->cd(2);
    hPIDgoodEventsWiderQave->DrawClone("colz");
    cPIDsPhysicalUnits->cd(3);
    hPIDnotGoodEventsQave->DrawClone("colz");
    cPIDsPhysicalUnits->cd(4);
    hPIDnotAtAllGoodEventsQave->DrawClone("colz");

    auto* cThetaPhi = new TCanvas("cThetaPhi", "Theta and Phi distributions", 1200, 800);
    cThetaPhi->Divide(2, 3);
    cThetaPhi->cd(1);
    hThetaGoodEvents->DrawClone();
    cThetaPhi->cd(2);
    hThetaAll->DrawClone();
    cThetaPhi->cd(3);
    hThetaBad->DrawClone();
    cThetaPhi->cd(4);
    hPhiGoodEvents->DrawClone();
    cThetaPhi->cd(5);
    hPhiAll->DrawClone();
    cThetaPhi->cd(6);
    hPhiBad->DrawClone();

    // Save bad events in a phi window, for cross-checking
    std::ofstream outFilePhi_130_135_BadEvents("./Outputs/eventsPhi_130_135_BadEvents.dat");
    defSummary.Foreach(
        [&outFilePhi_130_135_BadEvents, &cuts](ActRoot::MergerData& m, double maxGap, double meanGap)
        {
            if(m.fPhiLight > 130 && m.fPhiLight < 135 && !cuts.IsInside("goodEvents", maxGap, meanGap))
            {
                m.Stream(outFilePhi_130_135_BadEvents);
            }
        },
        {"MergerData", "maxGap", "meanGap"});
}