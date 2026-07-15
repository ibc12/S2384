// TracksMissingPadsXZYZProjection.C
//
// Studies missing-pad gaps in the XZ / YZ occupancy of the light cluster,
// including a "denser projection" choice per event and a restriction to the
// Bragg-peak region at the end of the track.
//
// This is one half of what used to be a single macro; the s-coordinate gap
// analysis (and the definition of the good/bad event quality cuts) now lives
// in TracksMissingPadsGapAnalysis.C. The maxGap/meanGap columns and the
// event-quality cuts are recomputed here too (they are cheap) so this macro
// can be run independently and still classify good/bad events the same way.
// Struct/function definitions are shared via GapAnalysisHelpers.h.

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

#include "MissingPads.h"

#include <iostream>
#include <string>

void TracksMissingPadsProjectionsAnalysis()
{
    ROOT::EnableImplicitMT(); // Multithreading is worthwhile since we process the full dataset.

    ActRoot::InputParser parser {};
    parser.ReadFile("../../configs/detector.conf");
    auto driftBlock = parser.GetBlock("Merger");
    auto driftFactor = driftBlock->GetDouble("DriftFactor");

    // Fraction of the track (counted backwards from the end, i.e. from max s)
    // used for the row-vs-Z gap-counting restricted to the Bragg-peak region.
    // Adjust this to whatever window makes physical sense; 0.5 keeps the
    // last 50% of the projected track length.
    double fracFromEndProj = 1;

    ROOT::RDataFrame df {"PreProcessed_Tree", "../../PostAnalysis/Outputs/tree_preprocess_F_11Li.root"};

    auto dfFiltered = df.Filter([](ActRoot::ModularData& m, ActRoot::MergerData& merger)
                                { return m.Get("GATCONF") == 8; }, {"ModularData", "MergerData"}); // only L1

    auto defSummary =
        dfFiltered
            // Recompute the s-projection gap summary (cheap) so this macro
            // can reuse the same maxGap/meanGap event-quality cuts as
            // TracksMissingPadsGapAnalysis.C without depending on it.
            .Define("projPointsMM",
                    [&driftFactor](ActRoot::MergerData& m, ActRoot::TPCData tpc)
                    {
                        std::vector<ProjectedPoint> empty;
                        auto lightIdx = m.fLightIdx;
                        if(lightIdx < 0 || lightIdx >= (int)tpc.fClusters.size())
                            return empty;
                        auto lightCl = tpc.fClusters[lightIdx];
                        auto rp = m.fRP;
                        return GetProjectedPositionsMM(lightCl, rp);
                    },
                    {"MergerData", "TPCData"})
            .Define("gapSummary", [](const std::vector<ProjectedPoint>& points) { return ComputeGapSummary(points); },
                    {"projPointsMM"})
            .Define("maxGap", [](const GapSummary& g) { return g.maxGap; }, {"gapSummary"})
            .Define("meanGap", [](const GapSummary& g) { return g.meanGap; }, {"gapSummary"})
            .Define("nVoxelsG", [](const GapSummary& g) { return g.nVoxels; }, {"gapSummary"})
            .Filter([](const GapSummary& g) { return g.valid; },
                    {"gapSummary"}) // Discard events with fewer than 3 voxels.
            // --- XZ / YZ pad-projection counts (global occupancy reference) -
            .Define("padSummary",
                    [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                    {
                        PadProjectionSummary empty;
                        auto lightIdx = m.fLightIdx;
                        if(lightIdx < 0 || lightIdx >= (int)tpc.fClusters.size())
                            return empty;
                        auto lightCl = tpc.fClusters[lightIdx];
                        auto v = lightCl.GetVoxels();
                        return ComputePadProjectionCounts(lightCl.GetVoxels(), 1);
                    },
                    {"MergerData", "TPCData"})
            .Define("nPadsXZ", [](const PadProjectionSummary& p) { return p.nPadsXZ; }, {"padSummary"})
            .Define("nPadsYZ", [](const PadProjectionSummary& p) { return p.nPadsYZ; }, {"padSummary"})
            .Define("nPadsXY", [](const PadProjectionSummary& p) { return p.nPadsXY; }, {"padSummary"})
            .Define("nPadsXZnorm", [](const PadProjectionSummary& p, const GapSummary& g)
                    { return (double)p.nPadsXZ / (double)g.sRange; }, {"padSummary", "gapSummary"})
            .Define("nPadsYZnorm", [](const PadProjectionSummary& p, const GapSummary& g)
                    { return (double)p.nPadsYZ / (double)g.sRange; }, {"padSummary", "gapSummary"})
            .Define("nPadsXYnorm", [](const PadProjectionSummary& p, const GapSummary& g)
                    { return (double)p.nPadsXY / (double)g.sRange; }, {"padSummary", "gapSummary"})
            .Define("ratioXZoverYZ", [](const PadProjectionSummary& p) { return p.ratioXZoverYZ; }, {"padSummary"})
            // --- gaps restricted to the Bragg-peak region, computed on -------
            // --- whichever projection (XZ or YZ) is wider for this event -----
            .Define("chosenGapResult",
                    [fracFromEndProj](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                    {
                        ChosenProjectionResult empty;

                        auto lightIdx = m.fLightIdx;

                        if(lightIdx < 0 || lightIdx >= (int)tpc.fClusters.size())
                            return empty;

                        const auto& lightCl = tpc.fClusters[lightIdx];

                        return ComputeGapsDenserProjection(lightCl.GetVoxels(), fracFromEndProj);
                    },
                    {"MergerData", "TPCData"})
            .Define("projGapResult", [](const ChosenProjectionResult& c) { return c.result; }, {"chosenGapResult"})
            .Define("usedXprojection", [](const ChosenProjectionResult& c) { return c.usedX; }, {"chosenGapResult"})
            .Define("nGapsProjection", [](const ProjGapResult& r) { return r.nGaps; }, {"projGapResult"})
            .Define("zCoverageFracProj", [](const ProjGapResult& r) { return r.zCoverageFrac; }, {"projGapResult"})
            .Define("nZSlicesCoveredProj", [](const ProjGapResult& r) { return r.nZSlicesCovered; }, {"projGapResult"})
            .Define("nZSlicesExpectedProj", [](const ProjGapResult& r) { return r.nZSlicesExpected; },
                    {"projGapResult"})
            .Define("meanPadsPerSliceProj", [](const ProjGapResult& r) { return r.meanPadsPerSlice; },
                    {"projGapResult"})
            .Define("medianPadsPerSliceProj", [](const ProjGapResult& r) { return r.medianPadsPerSlice; },
                    {"projGapResult"})
            .Define("stdPadsPerSliceProj", [](const ProjGapResult& r) { return r.stdPadsPerSlice; }, {"projGapResult"})
            .Define("nSinglePadSlicesProj", [](const ProjGapResult& r) { return r.nSinglePadSlices; },
                    {"projGapResult"})
            .Define("fracSinglePadSlicesProj", [](const ProjGapResult& r) { return r.fracSinglePadSlices; },
                    {"projGapResult"})
            .Filter([](const ProjGapResult& r) { return r.valid; }, {"projGapResult"})
            .Filter([](const PadProjectionSummary& p) { return p.valid; }, {"padSummary"});

    auto dfSummary1PadSlices =
        defSummary.Filter([](const ProjGapResult& r) { return r.nSinglePadSlices > 1; }, {"projGapResult"});

    // --- XZ / YZ pad-projection occupancy (global reference, both kept) ----
    auto hNPadsXZ =
        defSummary.Histo1D({"hNPadsXZ", "Occupied cells in XZ projection;nPadsXZ;Events", 200, 0, 400}, "nPadsXZ");
    auto hNPadsYZ =
        defSummary.Histo1D({"hNPadsYZ", "Occupied cells in YZ projection;nPadsYZ;Events", 200, 0, 400}, "nPadsYZ");

    // Direct XZ-vs-YZ comparison: if pads are lost symmetrically this should
    // sit close to the diagonal; systematic deviation from y=x points to one
    // plane/readout losing more pads in Z than the other.
    auto hNPadsXZvsYZ = defSummary.Histo2D(
        {"hNPadsXZvsYZ", "nPadsXZ vs nPadsYZ;nPadsYZ;nPadsXZ", 200, 0, 400, 200, 0, 400}, "nPadsYZ", "nPadsXZ");

    auto hRatioXZoverYZ =
        defSummary.Histo1D({"hRatioXZoverYZ", "nPadsXZ / nPadsYZ;ratio;Events", 200, 0, 5}, "ratioXZoverYZ");

    // Does the XZ/YZ occupancy correlate with the gap found along the track
    // direction (s-coordinate)? If large maxGap events also show a depleted
    // nPadsXZ or nPadsYZ (relative to nVoxelsG), that's strong evidence the
    // gap really is a missing-pads-in-Z effect rather than just a sparse but
    // healthy track.
    auto hGapVsNPadsXZ = defSummary.Histo2D(
        {"hGapVsNPadsXZ", "maxGap vs nPadsXZ;nPadsXZ;maxGap", 200, 0, 400, 100, 0, 10}, "nPadsXZ", "maxGap");
    auto hGapVsNPadsXZnorm = defSummary.Histo2D(
        {"hGapVsNPadsXZnorm", "maxGap vs nPadsXZ/nVoxelsG;nPadsXZ/nVoxelsG;maxGap", 200, 0, 1.5, 100, 0, 10},
        "nPadsXZnorm", "maxGap");
    auto hGapVsNPadsYZ = defSummary.Histo2D(
        {"hGapVsNPadsYZ", "maxGap vs nPadsYZ;nPadsYZ;maxGap", 200, 0, 400, 100, 0, 10}, "nPadsYZ", "maxGap");
    auto hGapVsNPadsYZnorm = defSummary.Histo2D(
        {"hGapVsNPadsYZnorm", "maxGap vs nPadsYZ/nVoxelsG;nPadsYZ/nVoxelsG;maxGap", 200, 0, 1.5, 100, 0, 10},
        "nPadsYZnorm", "maxGap");

    // Occupancy relative to raw voxel count: nPadsXZ (or YZ) / nVoxelsG close
    // to 1 means almost every voxel sits in its own XZ cell (sparse, likely
    // gappy); much less than 1 means many voxels share XZ cells (dense,
    // well-sampled track).
    auto defWithFillFrac =
        defSummary
            .Define("fillFracXZ", [](int nPadsXZ, int nVoxels)
                    { return nVoxels > 0 ? (double)nPadsXZ / nVoxels : -1.; }, {"nPadsXZ", "nVoxelsG"})
            .Define("fillFracYZ", [](int nPadsYZ, int nVoxels)
                    { return nVoxels > 0 ? (double)nPadsYZ / nVoxels : -1.; }, {"nPadsYZ", "nVoxelsG"});

    auto hFillFracXZ =
        defWithFillFrac.Histo1D({"hFillFracXZ", "nPadsXZ / nVoxelsG;fill fraction;Events", 100, 0, 1.5}, "fillFracXZ");
    auto hFillFracYZ =
        defWithFillFrac.Histo1D({"hFillFracYZ", "nPadsYZ / nVoxelsG;fill fraction;Events", 100, 0, 1.5}, "fillFracYZ");

    // --- which projection (XZ or YZ) was chosen for the gap analysis -------
    // usedXprojection == true -> XZ was wider (denser) and got used;
    // usedXprojection == false -> YZ was wider and got used.
    // This lets you check whether the choice is roughly 50/50 or strongly
    // tied to phi (as you'd expect if the track's "width" direction rotates
    // with the emission angle).
    auto hUsedProjection = defSummary.Histo1D(
        {"hUsedProjection", "Projection chosen for gap analysis;0 = YZ chosen, 1 = XZ chosen;Events", 2, -0.5, 1.5},
        "usedXprojection");
    auto hUsedProjectionVsPhi = defSummary.Histo2D(
        {"hUsedProjectionVsPhi", "Chosen projection vs phi;phi [deg];0 = YZ, 1 = XZ", 360, -180, 180, 2, -0.5, 1.5},
        "fPhiLight", "usedXprojection");
    auto hUsedProjectionVsTheta = defSummary.Histo2D(
        {"hUsedProjectionVsTheta", "Chosen projection vs theta;theta [deg];0 = YZ, 1 = XZ", 180, 0, 180, 2, -0.5, 1.5},
        "fThetaLight", "usedXprojection");

    // --- Pads per Z slice, on the chosen (denser) projection ----------------
    auto hMeanPadsPerSliceProj = defSummary.Histo1D(
        {"hMeanPadsPerSliceProj", "Mean pads per Z slice (chosen projection);Mean pads/slice;Events", 50, 0, 10},
        "meanPadsPerSliceProj");

    auto hMedianPadsPerSliceProj = defSummary.Histo1D(
        {"hMedianPadsPerSliceProj", "Median pads per Z slice (chosen projection);Median pads/slice;Events", 20, 0, 10},
        "medianPadsPerSliceProj");

    auto hStdPadsPerSliceProj =
        defSummary.Histo1D({"hStdPadsPerSliceProj",
                            "Std. dev. of pads per Z slice (chosen projection);#sigma(pads/slice);Events", 100, 0, 3},
                           "stdPadsPerSliceProj");

    auto hNSinglePadSlicesProj = defSummary.Histo1D(
        {"hNSinglePadSlicesProj", "Slices with one pad (chosen projection);N slices;Events", 20, 0, 20},
        "nSinglePadSlicesProj");

    auto hFracSinglePadSlicesProj = defSummary.Histo1D(
        {"hFracSinglePadSlicesProj", "Fraction of one-pad slices (chosen projection);Fraction;Events", 100, 0, 1},
        "fracSinglePadSlicesProj");

    auto hFracSingleVsMaxGap = defSummary.Histo2D(
        {"hFracSingleVsMaxGap", "Fraction of one-pad slices vs max gap;Max gap [mm];Fraction", 100, 0, 10, 100, 0, 1},
        "maxGap", "fracSinglePadSlicesProj");

    auto hFracSingleVsMeanGap =
        defSummary.Histo2D({"hFracSingleVsMeanGap", "Fraction of one-pad slices vs mean gap;Mean gap [mm];Fraction",
                            100, 0, 3.5, 100, 0, 1},
                           "meanGap", "fracSinglePadSlicesProj");

    auto hFracSingleVsNGapsProj =
        defSummary.Histo2D({"hFracSingleVsNGapsProj",
                            "Fraction of one-pad slices vs gaps (chosen proj.);N gaps;Fraction", 20, 0, 20, 100, 0, 1},
                           "nGapsProjection", "fracSinglePadSlicesProj");

    auto hMeanPadsVsMaxGap = defSummary.Histo2D(
        {"hMeanPadsVsMaxGap", "Mean pads/slice vs max gap;Max gap [mm];Mean pads/slice", 100, 0, 10, 60, 0, 6},
        "maxGap", "meanPadsPerSliceProj");

    auto hStdPadsVsMaxGap = defSummary.Histo2D(
        {"hStdPadsVsMaxGap", "Std pads/slice vs max gap;Max gap [mm];Std pads/slice", 100, 0, 10, 60, 0, 3}, "maxGap",
        "stdPadsPerSliceProj");

    auto hDeltaZ = dfFiltered
                       .Define("deltaZ",
                               [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                               {
                                   const auto& vox = tpc.fClusters[m.fLightIdx].GetVoxels();

                                   int zmin = 1e9;
                                   int zmax = -1e9;

                                   for(const auto& v : vox)
                                   {
                                       int z = std::lround(v.GetPosition().Z());
                                       zmin = std::min(zmin, z);
                                       zmax = std::max(zmax, z);
                                   }

                                   return zmax - zmin;
                               },
                               {"MergerData", "TPCData"})
                       .Histo1D({"hDeltaZ", "", 100, 0, 100}, "deltaZ");

    auto hMaxGapVSMeanGap1PadSlicesProj = dfSummary1PadSlices.Histo2D(
        {"hMaxGapVSMeanGap1PadSlicesProj", "maxGap vs meanGap (events with 1-pad slices);meanGap [mm];maxGap [mm]", 100,
         0, 10, 100, 0, 3.5},
        "maxGap", "meanGap");
    auto hPhi1PadSlicesProj = dfSummary1PadSlices.Histo1D(
        {"hPhi1PadSlicesProj", "Phi of events with 1-pad slices;phi [deg];Events", 360, -180, 180}, "fPhiLight");
    auto hStdPads1PadSlicesProj = dfSummary1PadSlices.Histo1D(
        {"hStdPads1PadSlicesProj", "Std. dev. of pads/slice (events with 1-pad slices);#sigma(pads/slice);Events", 100,
         0, 10},
        "stdPadsPerSliceProj");

    // --- Canvases for XZ / YZ pad projections (global occupancy) ----------
    auto* c4 = new TCanvas("c4", "XZ / YZ pad occupancy", 1200, 800);
    c4->Divide(2, 2);
    c4->cd(1);
    hNPadsXZ->DrawClone();
    c4->cd(2);
    hNPadsYZ->DrawClone();
    c4->cd(3);
    hNPadsXZvsYZ->DrawClone("colz");
    c4->cd(4);
    hRatioXZoverYZ->DrawClone();

    auto* c5 = new TCanvas("c5", "Max Gap vs XZ/YZ occupancy", 1200, 400);
    c5->Divide(2, 2);
    c5->cd(1);
    hGapVsNPadsXZ->DrawClone("colz");
    c5->cd(2);
    hGapVsNPadsYZ->DrawClone("colz");
    c5->cd(3);
    hGapVsNPadsXZnorm->DrawClone("colz");
    c5->cd(4);
    hGapVsNPadsYZnorm->DrawClone("colz");

    auto* c6 = new TCanvas("c6", "XZ / YZ fill fraction", 800, 400);
    c6->Divide(2, 1);
    c6->cd(1);
    hFillFracXZ->DrawClone();
    c6->cd(2);
    hFillFracYZ->DrawClone();

    // --- which projection got picked for the gap analysis -----------------
    auto* cUsedProjection = new TCanvas("cUsedProjection", "Chosen projection (XZ vs YZ) for gap analysis", 1500, 500);
    cUsedProjection->Divide(3, 1);
    cUsedProjection->cd(1);
    hUsedProjection->DrawClone();
    cUsedProjection->cd(2);
    hUsedProjectionVsPhi->DrawClone("colz");
    cUsedProjection->cd(3);
    hUsedProjectionVsTheta->DrawClone("colz");

    auto* cPadsPerSlice = new TCanvas("cPadsPerSlice", "Pads per Z slice (chosen projection)", 1500, 800);
    cPadsPerSlice->Divide(3, 2);
    cPadsPerSlice->cd(1);
    hMeanPadsPerSliceProj->DrawClone();
    cPadsPerSlice->cd(2);
    hMedianPadsPerSliceProj->DrawClone();
    cPadsPerSlice->cd(3);
    hStdPadsPerSliceProj->DrawClone();
    cPadsPerSlice->cd(4);
    hNSinglePadSlicesProj->DrawClone();
    cPadsPerSlice->cd(5);
    hFracSinglePadSlicesProj->DrawClone();

    auto* cPadsPerSliceCorr =
        new TCanvas("cPadsPerSliceCorr", "Pads per slice correlations (chosen projection)", 1500, 800);
    cPadsPerSliceCorr->Divide(3, 2);
    cPadsPerSliceCorr->cd(1);
    hFracSingleVsMaxGap->DrawClone("colz");
    cPadsPerSliceCorr->cd(2);
    hFracSingleVsMeanGap->DrawClone("colz");
    cPadsPerSliceCorr->cd(3);
    hFracSingleVsNGapsProj->DrawClone("colz");
    cPadsPerSliceCorr->cd(4);
    hMeanPadsVsMaxGap->DrawClone("colz");
    cPadsPerSliceCorr->cd(5);
    hStdPadsVsMaxGap->DrawClone("colz");

    auto* cDeltaZ = new TCanvas("cDeltaZ", "Z extent of the light cluster", 800, 600);
    hDeltaZ->DrawClone();

    auto* c1PadSlicesProj = new TCanvas("c1PadSlicesProj", "Events with 1-pad slices (chosen projection)", 1200, 400);
    c1PadSlicesProj->Divide(3, 1);
    c1PadSlicesProj->cd(1);
    hMaxGapVSMeanGap1PadSlicesProj->DrawClone("colz");
    c1PadSlicesProj->cd(2);
    hPhi1PadSlicesProj->DrawClone();
    c1PadSlicesProj->cd(3);
    hStdPads1PadSlicesProj->DrawClone();

    std::cout << "Total events processed: " << *defSummary.Count() << std::endl;

    // Reuse the same event-quality cuts as TracksMissingPadsGapAnalysis.C
    // (defined on maxGap vs meanGap) to classify good/bad events for the
    // nGapsProjection / zCoverage plots below.
    ActRoot::CutsManager<std::string> cuts {};
    cuts.ReadCut("goodEvents", "./Cuts/cut_goodEvents_maxGap_meanGap.root");

    auto dfGoodEvents =
        defSummary.Filter([&cuts](double meanGap, double maxGap)
                          { return cuts.IsInside("goodEvents", maxGap, meanGap); }, {"meanGap", "maxGap"});

    // Plots for the nGaps in the chosen (denser) projection (Bragg-peak region only)
    auto hNGapsProjGoodEvents = dfGoodEvents.Histo1D(
        {"hNGapsProjGoodEvents", "Number of gaps (chosen projection, Bragg region) for good events;nGaps;Events", 20, 0,
         20},
        "nGapsProjection");
    auto hNGapsProjAll = defSummary.Histo1D(
        {"hNGapsProjAll", "Number of gaps (chosen projection, Bragg region) for all events;nGaps;Events", 20, 0, 20},
        "nGapsProjection");
    auto hMaxGapVSMeanGapNGapsProj =
        defSummary.Filter([](int& nGapsProj) { return nGapsProj > 0; }, {"nGapsProjection"})
            .Histo2D({"hMaxGapVSMeanGapNGapsProj", "Max gap vs mean gap for nGaps > 0;max gap [mm];mean gap [mm]", 200,
                      0, 10, 200, 0, 3.5},
                     "maxGap", "meanGap");
    auto hMaxGapVSMeanNoGapNGapsProj =
        defSummary.Filter([](int& nGapsProj) { return nGapsProj == 0; }, {"nGapsProjection"})
            .Histo2D({"hMaxGapVSMeanNoGapNGapsProj", "Max gap vs mean gap for nGaps == 0;max gap [mm];mean gap [mm]",
                      200, 0, 10, 200, 0, 3.5},
                     "maxGap", "meanGap");

    // --- z-coverage of the Bragg-peak window, on the chosen projection -----
    // Tells you what fraction of the z-slices spanned by the last
    // fracFromEndProj of the track are actually populated. A low value plus a
    // high nGapsProjection points to a genuinely depleted region, not just
    // one or two isolated missing pads.
    auto hZCoverageFracProj =
        defSummary.Histo1D({"hZCoverageFracProj",
                            "Z-slice coverage in Bragg-peak region (chosen proj.);covered/expected;Events", 100, 0, 1},
                           "zCoverageFracProj");
    auto hGapsVsCoverageProj = defSummary.Histo2D(
        {"hGapsVsCoverageProj", "nGaps vs z-coverage (chosen proj.);covered/expected;nGaps", 100, 0, 1, 20, 0, 20},
        "zCoverageFracProj", "nGapsProjection");

    auto* cNGapsProj = new TCanvas("cNGapsProj", "Number of gaps (chosen projection, Bragg region)", 800, 600);
    cNGapsProj->Divide(2, 2);
    cNGapsProj->cd(1);
    hNGapsProjGoodEvents->DrawClone();
    cNGapsProj->cd(2);
    hNGapsProjAll->DrawClone();
    cNGapsProj->cd(3);
    hMaxGapVSMeanGapNGapsProj->DrawClone("colz");
    cNGapsProj->cd(4);
    hMaxGapVSMeanNoGapNGapsProj->DrawClone("colz");

    auto* cZCoverageProj = new TCanvas("cZCoverageProj", "Z coverage in Bragg-peak region (chosen proj.)", 800, 400);
    cZCoverageProj->Divide(2, 1);
    cZCoverageProj->cd(1);
    hZCoverageFracProj->DrawClone();
    cZCoverageProj->cd(2);
    hGapsVsCoverageProj->DrawClone("colz");

    // Save some events if needed
    // std::ofstream outFileNGapsProj_NoGapsBad("./Outputs/eventsNGapsProj_NoGapsBad.dat");
    // std::ofstream outFileNGapsProj_GapsGood("./Outputs/eventsNGapsProj_GapsGood.dat");
    // std::ofstream outFileNGapsProj_Events1PadSlicesProj("./Outputs/eventsNGapsProj_Events1PadSlicesProj.dat");
    // defSummary.Foreach(
    //     [&outFileNGapsProj_Events1PadSlicesProj](ActRoot::MergerData& m, double fracSinglePadSlices)
    //     {
    //         if(fracSinglePadSlices > 0)
    //         {
    //             m.Stream(outFileNGapsProj_Events1PadSlicesProj);
    //         }
    //     },
    //     {"MergerData", "fracSinglePadSlicesProj"});
}