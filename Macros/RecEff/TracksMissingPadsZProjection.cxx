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

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <set>
#include <string>
#include <utility>
#include <vector>


struct ProjectedPoint
{
    double s;
    ROOT::Math::XYZPointF pos3D;
};

// Same as before: project the voxels onto the cluster line and sort them by s.
std::vector<ProjectedPoint>
GetProjectedPositions(ActRoot::Cluster lightCl, const ROOT::Math::XYZPointF& ref, double driftFactor = 1.)
{
    // auto voxels = lightCl.GetRefToVoxels();
    auto line = lightCl.GetLine();
    line.Scale(2, driftFactor);
    lightCl.ScaleVoxels(2, driftFactor);
    auto voxels = lightCl.GetRefToVoxels();
    auto dir = line.GetDirection().Unit();
    lightCl.SortAlongDir(dir);

    std::vector<ProjectedPoint> points;
    points.reserve(voxels.size());
    for(auto& v : lightCl.GetVoxels())
    {
        auto pos = v.GetPosition();
        auto projPoint = line.ProjectionPointOnLine(pos);
        auto diff = projPoint - ref;
        double s = diff.X() * dir.X() + diff.Y() * dir.Y() + diff.Z() * dir.Z();
        points.push_back({s, pos});
    }
    std::sort(points.begin(), points.end(), [](const ProjectedPoint& a, const ProjectedPoint& b) { return a.s < b.s; });
    return points;
}

// Summary of the gaps in a single event: maximum gap (absolute and normalized),
// number of voxels, and total track length in s.
//
// The normalization is important: a 5 mm gap can be huge in a short, dense
// track but completely negligible in a long, sparse one. Therefore, we divide
// by the "typical" gap (the median of the event's own gaps) to obtain a
// dimensionless quantity that can be compared across events.
struct GapSummary
{
    double maxGap = -1;
    double maxGapNorm = -1; // maxGap / median(gaps)
    double medianGap = -1;
    double meanGap = -1;
    double stdGap = -1;
    double sRange = -1;
    int nVoxels = 0;
    bool valid = false;
};

GapSummary ComputeGapSummary(const std::vector<ProjectedPoint>& points)
{
    GapSummary res;
    res.nVoxels = (int)points.size();
    if(points.size() < 3) // With fewer than 3 points there is no reliable "typical" gap.
        return res;

    std::vector<double> gaps;
    gaps.reserve(points.size() - 1);
    for(size_t i = 1; i < points.size(); ++i)
        gaps.push_back(points[i].s - points[i - 1].s);

    auto maxIt = std::max_element(gaps.begin(), gaps.end());
    res.maxGap = *maxIt;

    std::vector<double> sortedGaps = gaps;
    std::sort(sortedGaps.begin(), sortedGaps.end());
    res.medianGap = sortedGaps[sortedGaps.size() / 2];
    res.meanGap = std::accumulate(sortedGaps.begin(), sortedGaps.end(), 0.0) / sortedGaps.size();

    // if(res.medianGap > 0)
    res.maxGapNorm = res.maxGap / res.medianGap;

    // Calculate standard deviation of gaps
    double sumSquaredDiffs = 0;
    for(const auto& gap : gaps)
    {
        double diff = gap - res.meanGap;
        sumSquaredDiffs += diff * diff;
    }
    res.stdGap = std::sqrt(sumSquaredDiffs / gaps.size());

    res.sRange = points.back().s - points.front().s;
    res.valid = true;
    return res;
}

// --- New: pad counts in the XZ and YZ projections -------------------------
//
// Idea: a voxel is really "pad (X,Y) fired at drift time Z". If pads are being
// lost specifically in Z (e.g. dead time buckets, saturation, thresholding in
// the drift coordinate) that loss should show up as gaps when you look at the
// (X,Z) and (Y,Z) occupancy separately, and it may be asymmetric between the
// two projections if one pad plane/readout is behaving differently from the
// other. Counting the number of *distinct occupied cells* (rather than raw
// voxel counts) tells you how much of the XZ / YZ footprint is actually
// filled versus how sparse/gappy it is once you collapse the third dimension.
struct PadProjectionSummary
{
    int nPadsXZ = 0;
    int nPadsYZ = 0;
    int nPadsXY = 0; // occupancy in the "real" pad plane, useful as a reference
    double ratioXZoverYZ = -1;
    bool valid = false;
};

// padPitch: bin size (mm) used to group nearby positions into the same "pad
// cell" in the projection. Set it to your actual pad pitch (commonly 2 mm for
// ACTAR TPC); if voxel positions are already discretized to pad centers, any
// pitch <= that spacing works fine as a binning tolerance.
PadProjectionSummary ComputePadProjectionCounts(const std::vector<ActRoot::Voxel>& voxels, double padPitch = 2.0)
{
    PadProjectionSummary res;
    if(voxels.empty() || padPitch <= 0)
        return res;

    std::set<std::pair<int, int>> xz;
    std::set<std::pair<int, int>> yz;
    std::set<std::pair<int, int>> xy;

    for(const auto& v : voxels)
    {
        auto pos = v.GetPosition();
        int ix = (int)std::round(pos.X() / padPitch);
        int iy = (int)std::round(pos.Y() / padPitch);
        int iz = (int)std::round(pos.Z() / padPitch);

        xz.insert({ix, iz});
        yz.insert({iy, iz});
        xy.insert({ix, iy});
    }

    res.nPadsXZ = (int)xz.size();
    res.nPadsYZ = (int)yz.size();
    res.nPadsXY = (int)xy.size();
    if(res.nPadsYZ > 0)
        res.ratioXZoverYZ = (double)res.nPadsXZ / (double)res.nPadsYZ;
    res.valid = true;
    return res;
}

void TracksMissingPadsZProjection()
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

    // Pad pitch used for the XZ/YZ/XY occupancy counting below.

    // ROOT::RDataFrame df {*chain};

    ROOT::RDataFrame df {"PreProcessed_Tree", "../../PostAnalysis/Outputs/tree_preprocess_F_11Li.root"};

    auto dfFiltered = df.Filter([](ActRoot::ModularData& m, ActRoot::MergerData& merger)
                                { return m.Get("GATCONF") == 8; }, {"ModularData", "MergerData"}); // only L1

    // Previously done with && merger.fRun == 69

    // No event selection: process the entire chain.
    // If needed, any physics selection (e.g. only L1 events) can be added here.
    auto defSummary =
        dfFiltered
            .Define("gapSummary",
                    [&driftFactor](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                    {
                        GapSummary empty;
                        auto lightIdx = m.fLightIdx;
                        if(lightIdx < 0 || lightIdx >= (int)tpc.fClusters.size())
                            return empty;
                        auto lightCl = tpc.fClusters[lightIdx];
                        auto rp = m.fRP;
                        auto points = GetProjectedPositions(lightCl, rp, driftFactor);
                        return ComputeGapSummary(points);
                    },
                    {"MergerData", "TPCData"})
            .Define("maxGap", [](const GapSummary& g) { return g.maxGap; }, {"gapSummary"})
            .Define("maxGapNorm", [](const GapSummary& g) { return g.maxGapNorm; }, {"gapSummary"})
            .Define("medianGap", [](const GapSummary& g) { return g.medianGap; }, {"gapSummary"})
            .Define("meanGap", [](const GapSummary& g) { return g.meanGap; }, {"gapSummary"})
            .Define("nVoxelsG", [](const GapSummary& g) { return g.nVoxels; }, {"gapSummary"})
            .Define("sRange", [](const GapSummary& g) { return g.sRange; }, {"gapSummary"})
            .Define("stdGap", [](const GapSummary& g) { return g.stdGap; }, {"gapSummary"})
            .Filter([](const GapSummary& g) { return g.valid; },
                    {"gapSummary"}) // Discard events with fewer than 3 voxels.
            // --- New: XZ / YZ pad-projection counts -----------------------
            .Define("padSummary",
                    [&driftFactor](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                    {
                        PadProjectionSummary empty;
                        auto lightIdx = m.fLightIdx;
                        if(lightIdx < 0 || lightIdx >= (int)tpc.fClusters.size())
                            return empty;
                        auto lightCl = tpc.fClusters[lightIdx];
                        // Apply the same Z scaling as GetProjectedPositions so the
                        // pad pitch is meaningful in the same units.
                        // lightCl.ScaleVoxels(2, driftFactor);
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
            .Filter([](const PadProjectionSummary& p) { return p.valid; }, {"padSummary"});

    // Get the data in the two branches of phi
    auto defSummaryPhiPositive =
        defSummary.Filter([](ActRoot::MergerData& m) { return m.fPhiLight > 0; }, {"MergerData"});
    auto defSummaryPhiNegative =
        defSummary.Filter([](ActRoot::MergerData& m) { return m.fPhiLight < 0; }, {"MergerData"});

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


    // --- Distribution of the normalized maximum gap (dimensionless) ---
    // This is likely the most useful observable for defining a threshold.
    // A maxGapNorm around 1–2 corresponds to the expected spacing, whereas
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

    // --- New: XZ / YZ pad-projection occupancy ---------------------------
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

    // Does the XZ/YZ occupancy correlate with the gap you already found along
    // the track direction? If large maxGapNorm events also show a depleted
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

    auto* c1 = new TCanvas("c1", "Max gap distributions", 1200, 800);
    c1->Divide(2, 2);
    c1->cd(1);
    // Fit to two gausians to see if there is a clear separation between the "normal" and "suspicious" populations.
    // auto gaus1 = new TF1("gaus1", "gaus", 0, 4);
    // auto gaus2 = new TF1("gaus2", "gaus", 0, 4);
    // auto gaus3 = new TF1("gaus", "gaus", 4, 10);
    // // Give initial parameters to the fits based on the histogram content
    // gaus1->SetParameters(150, 1.7, 0.5); // amplitude, mean, sigma
    // gaus2->SetParameters(400, 3, 0.5);
    // gaus3->SetParameters(20, 5.6, 0.5);
    // // Put bounds to the parameters
    // gaus1->SetParLimits(0, 100, 1e6); // amplitude between 0 and 1e6
    // gaus2->SetParLimits(0, 200, 1e6);
    // gaus1->SetParLimits(1, 1.5, 2); // mean between 0 and 10
    // gaus2->SetParLimits(1, 2.6, 3.4); // mean between 10 and 50
    // gaus3->SetParLimits(1, 5.3, 6);   // mean between 0 and 50
    // gaus1->SetParLimits(2, 0.1, 1); // sigma between 0 and 5
    // gaus2->SetParLimits(2, 0.1, 1);
    // gaus3->SetParLimits(2, 0.1, 1);

    // gaus1->FixParameter(1, 2.4);
    // gaus2->FixParameter(1, 3);
    // Fit the histogram with the three gaussians
    // hMaxGap->Fit(gaus1, "R");
    // hMaxGap->Fit(gaus2, "R+");
    // hMaxGap->Fit(gaus3, "R+");
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
    c33->Divide(2, 2);
    c33->cd(1);
    hMeanGap->DrawClone();
    c33->cd(2);
    hMeanGapVSMaxGap->DrawClone("colz");
    c33->cd(3);
    hMeanGapVSMaxGapPhiPositive->DrawClone("colz");
    c33->cd(4);
    hMeanGapVSMaxGapPhiNegative->DrawClone("colz");

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

    // --- New canvases for XZ / YZ pad projections -------------------------
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
            .Histo2D({"h    PIDallPhiPositive", "PID for all events (phi>0);TL;Qtot", 200, 0, 120, 2000, 0, 3e5},
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
}