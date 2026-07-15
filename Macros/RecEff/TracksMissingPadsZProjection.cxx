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
#include <functional>
#include <iostream>
#include <map>
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
    double mad = -1;      // median absolute deviation of gaps
    double maxGapZ = -1;  // (maxGap - medianGap) / (1.4826 * mad) -> robust z-score
    double skewness = -1; // third standardized moment of the gaps
    double sRange = -1;
    int nVoxels = 0;
    bool valid = true;
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

    // --- Robust z-score of maxGap using the MAD ---------------------------
    // MAD barely moves even if maxGap itself is a huge outlier, unlike stdGap.
    std::vector<double> absDev;
    absDev.reserve(gaps.size());
    for(const auto& gap : gaps)
        absDev.push_back(std::abs(gap - res.medianGap));
    std::sort(absDev.begin(), absDev.end());
    res.mad = absDev[absDev.size() / 2];

    double robustSigma = 1.4826 * res.mad; // consistency factor for a Gaussian
    res.maxGapZ = (robustSigma > 0) ? (res.maxGap - res.medianGap) / robustSigma : -1;

    // --- Skewness of the gap distribution ----------------------------------
    // Third moment: much more sensitive to a single extreme gap than stdGap.
    double sumCubedDiffs = 0;
    for(const auto& gap : gaps)
    {
        double diff = gap - res.meanGap;
        sumCubedDiffs += diff * diff * diff;
    }
    double m3 = sumCubedDiffs / gaps.size();
    res.skewness = (res.stdGap > 0) ? m3 / (res.stdGap * res.stdGap * res.stdGap) : -1;

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
    bool valid = true;
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

// --- Gaps in a "row" coordinate per Z-slice, restricted to a window at the
// END of the track (largest s), where the Bragg peak deposits most of the
// charge. -------------------------------------------------------------------
//
// This is now GENERIC: it doesn't hardcode X or Y. You pass in a lambda
// (getRow) that extracts whichever coordinate you want to use as the "row"
// grouped by Z (either pos.X() for the XZ projection or pos.Y() for the YZ
// projection). For each fixed z (a "row" of pads at a given drift time)
// *within that window*, take the sorted, distinct row-coordinate values that
// are occupied. Any place where two consecutive occupied values are not
// adjacent (diff > 1 pad unit) counts as ONE hole, no matter how many pads
// are missing in between.
// Example: at z=120, occupied row = {63, 64, 67} -> exactly one hole, since
// 65 and 66 are both missing between the 64 and 67 hits.
// Summing these holes over every z slice inside the window gives nGaps.
//
// We also report how much of the z-range in that window is actually covered
// by at least one hit (zCoverageFrac), since a handful of gaps in a densely
// sampled window means something very different than the same number of
// gaps in a window that's mostly empty in z.
struct ProjGapResult
{
    int nGaps = 0;

    int nZSlicesCovered = 0;
    int nZSlicesExpected = 0;
    double zCoverageFrac = 0.;

    double meanPadsPerSlice = 0.;
    double medianPadsPerSlice = 0.;
    double stdPadsPerSlice = 0.;

    int nSinglePadSlices = 0;
    double fracSinglePadSlices = 0.;

    bool valid = true;
};

ProjGapResult ComputeGapsProjectionInRange(const std::vector<ActRoot::Voxel>& voxels,
                                            const std::function<int(const ROOT::Math::XYZPointF&)>& getRow,
                                            double fracFromEnd = 1., int minDeltaZ = 0, int maxHoleSize = 2)
{
    ProjGapResult res;

    if(voxels.size() < 2)
        return res;

    if(fracFromEnd <= 0. || fracFromEnd > 1.)
        return res;

    //------------------------------------------------------------
    // Determine Z range
    //------------------------------------------------------------

    int zMin = std::numeric_limits<int>::max();
    int zMax = std::numeric_limits<int>::lowest();

    for(const auto& voxel : voxels)
    {
        auto pos = voxel.GetPosition();

        int iz = std::lround(pos.Z());

        zMin = std::min(zMin, iz);
        zMax = std::max(zMax, iz);
    }

    int deltaZ = zMax - zMin;

    if(deltaZ < minDeltaZ)
        return res;

    //------------------------------------------------------------
    // Keep only the last fracFromEnd of the Z range
    //------------------------------------------------------------

    int zThreshold = std::lround(zMax - fracFromEnd * deltaZ);

    std::map<int, std::set<int>> zToRow;

    for(const auto& voxel : voxels)
    {
        auto pos = voxel.GetPosition();

        int iz = std::lround(pos.Z());

        if(iz < zThreshold)
            continue;

        zToRow[iz].insert(getRow(pos));
    }

    if(zToRow.empty())
        return res;

    //------------------------------------------------------------
    // Count row gaps
    //------------------------------------------------------------

    std::vector<int> occupancies;
    occupancies.reserve(zToRow.size());
    int nGaps = 0;

    for(const auto& [z, rows] : zToRow)
    {
        occupancies.push_back((int)rows.size());

        if(rows.size() == 1)
            ++res.nSinglePadSlices;

        if(rows.size() < 2)
            continue;

        auto it = rows.begin();
        int previous = *it;
        ++it;

        for(; it != rows.end(); ++it)
        {
            int holeSize = *it - previous - 1;

            if(holeSize >= 1 && holeSize <= maxHoleSize)
                ++nGaps;

            previous = *it;
        }
    }

    //------------------------------------------------------------
    // Z coverage
    //------------------------------------------------------------

    int firstZ = zToRow.begin()->first;
    int lastZ = zToRow.rbegin()->first;

    res.nZSlicesExpected = lastZ - firstZ + 1;
    res.nZSlicesCovered = (int)zToRow.size();

    if(res.nZSlicesExpected > 0)
        res.zCoverageFrac = (double)res.nZSlicesCovered / res.nZSlicesExpected;

    res.nGaps = nGaps;
    res.valid = true;

    //------------------------------------------------------------
    // Occupancy statistics
    //------------------------------------------------------------

    if(!occupancies.empty())
    {
        // Mean
        double sum = std::accumulate(occupancies.begin(), occupancies.end(), 0.0);
        res.meanPadsPerSlice = sum / occupancies.size();

        // Median
        auto occSorted = occupancies;
        std::sort(occSorted.begin(), occSorted.end());

        if(occSorted.size() % 2 == 0)
            res.medianPadsPerSlice = 0.5 * (occSorted[occSorted.size() / 2 - 1] + occSorted[occSorted.size() / 2]);
        else
            res.medianPadsPerSlice = occSorted[occSorted.size() / 2];

        // Standard deviation
        double var = 0.0;

        for(auto n : occupancies)
            var += (n - res.meanPadsPerSlice) * (n - res.meanPadsPerSlice);

        res.stdPadsPerSlice = std::sqrt(var / occupancies.size());

        // Fraction of slices with only one pad
        res.fracSinglePadSlices = static_cast<double>(res.nSinglePadSlices) / occupancies.size();
    }

    return res;
}

// --- New: choose, per event, whichever projection (XZ or YZ) is "wide" ----
//
// A track that is thin in YZ but wide in XZ (or vice versa) should NOT
// always be analyzed in XZ: fixing the projection a priori mixes up "real
// missing pads in Z" with "the track is just narrow in this projection by
// geometry". Instead, compute the gap analysis in BOTH projections and keep
// the one with the larger mean pads/slice (i.e. the direction in which the
// track actually has width to lose pads from). usedX tells you which one
// was picked, so you can cross-check it against e.g. phi.
struct ChosenProjectionResult
{
    ProjGapResult result;
    bool usedX = true; // true -> XZ projection chosen, false -> YZ projection chosen
};

ChosenProjectionResult ComputeGapsDenserProjection(const std::vector<ActRoot::Voxel>& voxels, double fracFromEnd = 0.5,
                                                    int minDeltaZ = 30, int maxHoleSize = 2)
{
    auto getX = [](const ROOT::Math::XYZPointF& pos) { return (int)std::lround(pos.X()); };
    auto getY = [](const ROOT::Math::XYZPointF& pos) { return (int)std::lround(pos.Y()); };

    auto resX = ComputeGapsProjectionInRange(voxels, getX, fracFromEnd, minDeltaZ, maxHoleSize);
    auto resY = ComputeGapsProjectionInRange(voxels, getY, fracFromEnd, minDeltaZ, maxHoleSize);

    ChosenProjectionResult chosen;

    if(!resX.valid && !resY.valid)
    {
        chosen.result = resX; // invalid, valid == false
        return chosen;
    }
    if(resX.valid && !resY.valid)
    {
        chosen.result = resX;
        chosen.usedX = true;
        return chosen;
    }
    if(!resX.valid && resY.valid)
    {
        chosen.result = resY;
        chosen.usedX = false;
        return chosen;
    }

    // Both valid: pick the projection with the larger mean pads/slice, i.e.
    // the "wide" direction, since that's the one where a real gap in Z shows
    // up cleanly instead of being an artifact of the track just being thin
    // in that projection.
    if(resX.meanPadsPerSlice >= resY.meanPadsPerSlice)
    {
        chosen.result = resX;
        chosen.usedX = true;
    }
    else
    {
        chosen.result = resY;
        chosen.usedX = false;
    }
    return chosen;
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

    // Fraction of the track (counted backwards from the end, i.e. from max s)
    // used for the row-vs-Z gap-counting restricted to the Bragg-peak region.
    // Adjust this to whatever window makes physical sense; 0.3 keeps the
    // last 30% of the projected track length.
    double fracFromEndProj = 0.5;

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
            // Compute the projected/sorted points ONCE and reuse them both for
            // gapSummary and for the row-vs-Z gap-in-Bragg-peak-region calculation.
            .Define("projPoints",
                    [&driftFactor](ActRoot::MergerData& m, ActRoot::TPCData tpc)
                    {
                        std::vector<ProjectedPoint> empty;
                        auto lightIdx = m.fLightIdx;
                        if(lightIdx < 0 || lightIdx >= (int)tpc.fClusters.size())
                            return empty;
                        auto lightCl = tpc.fClusters[lightIdx];
                        auto rp = m.fRP;
                        return GetProjectedPositions(lightCl, rp, driftFactor);
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
                    {"gapSummary"}) // Discard events with fewer than 3 voxels.
            // --- XZ / YZ pad-projection counts (global occupancy reference) -
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
            // --- New: gaps restricted to the Bragg-peak region, computed on --
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

    auto minFrac = *defSummary.Min("fracSinglePadSlicesProj");

    auto maxFrac = *defSummary.Max("fracSinglePadSlicesProj");

    std::cout << minFrac << " " << maxFrac << std::endl;

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

    // --- New: which projection (XZ or YZ) was chosen for the gap analysis --
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

    auto hStdPadsPerSliceProj = defSummary.Histo1D(
        {"hStdPadsPerSliceProj", "Std. dev. of pads per Z slice (chosen projection);#sigma(pads/slice);Events", 100, 0,
         3},
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

    auto hFracSingleVsNGapsProj = defSummary.Histo2D(
        {"hFracSingleVsNGapsProj", "Fraction of one-pad slices vs gaps (chosen proj.);N gaps;Fraction", 20, 0, 20,
         100, 0, 1},
        "nGapsProjection", "fracSinglePadSlicesProj");

    auto hMeanPadsVsMaxGap = defSummary.Histo2D(
        {"hMeanPadsVsMaxGap", "Mean pads/slice vs max gap;Max gap [mm];Mean pads/slice", 100, 0, 10, 60, 0, 6},
        "maxGap", "meanPadsPerSliceProj");

    auto hStdPadsVsMaxGap = defSummary.Histo2D(
        {"hStdPadsVsMaxGap", "Std pads/slice vs max gap;Max gap [mm];Std pads/slice", 100, 0, 10, 60, 0, 3}, "maxGap",
        "stdPadsPerSliceProj");

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

    // --- New canvas: which projection got picked for the gap analysis -----
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

    auto* cPadsPerSliceCorr = new TCanvas("cPadsPerSliceCorr", "Pads per slice correlations (chosen projection)", 1500,
                                          800);
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

    // Plots for the nGaps in the chosen (denser) projection (Bragg-peak region only)
    auto hNGapsProjGoodEvents = dfGoodEvents.Histo1D(
        {"hNGapsProjGoodEvents", "Number of gaps (chosen projection, Bragg region) for good events;nGaps;Events", 20,
         0, 20},
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
    auto hZCoverageFracProj = defSummary.Histo1D(
        {"hZCoverageFracProj", "Z-slice coverage in Bragg-peak region (chosen proj.);covered/expected;Events", 100, 0,
         1},
        "zCoverageFracProj");
    auto hGapsVsCoverageProj = defSummary.Histo2D(
        {"hGapsVsCoverageProj", "nGaps vs z-coverage (chosen proj.);covered/expected;nGaps", 100, 0, 1, 20, 0, 20},
        "zCoverageFracProj", "nGapsProjection");

    // Theta and phi dstribution for good bad and all events
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

    // Save some events if needed
    // std::ofstream outFileNGapsProj_NoGapsBad("./Outputs/eventsNGapsProj_NoGapsBad.dat");
    // std::ofstream outFileNGapsProj_GapsGood("./Outputs/eventsNGapsProj_GapsGood.dat");
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