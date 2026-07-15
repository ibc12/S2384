#ifndef MissingPads_h
#define MissingPads_h

// Shared helpers for the two gap-analysis macros:
//   - TracksMissingPadsGapAnalysis.C   (gaps along the projected s-coordinate)
//   - TracksMissingPadsXZYZProjection.C (gaps in the XZ / YZ pad projections)
//
// Keeping these structs/functions in one header avoids duplicating the
// analysis logic while still letting each macro recompute the (cheap)
// quantities it needs from the other domain (e.g. macro 2 recomputing
// maxGap/meanGap to reuse the same event-quality cuts).

#include "ActCluster.h"
#include "ActLine.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActVoxel.h"

#include "Math/Point3Dfwd.h"
#include "Math/Vector3D.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <utility>
#include <vector>

// ============================================================================
// Section 1: projection onto the cluster line + gap statistics along s
// ============================================================================

struct ProjectedPoint
{
    double s;
    ROOT::Math::XYZPointF pos3D;
};

// Project the voxels onto the cluster line and sort them by s.
inline std::vector<ProjectedPoint>
GetProjectedPositionsMM(ActRoot::Cluster lightCl, const ROOT::Math::XYZPointF& ref, double driftFactor = 1.)
{
    auto line = lightCl.GetLine();
    line.Scale(2, driftFactor);
    lightCl.ScaleVoxels(2, driftFactor);
    auto dir = line.GetDirection().Unit();
    lightCl.SortAlongDir(dir);

    std::vector<ProjectedPoint> points;
    points.reserve(lightCl.GetVoxels().size());
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

inline std::vector<ProjectedPoint> GetProjectedPositionsPads(ActRoot::Cluster lightCl, const ROOT::Math::XYZPointF& ref)
{
    auto line = lightCl.GetLine();
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

inline GapSummary ComputeGapSummary(const std::vector<ProjectedPoint>& points)
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

// ============================================================================
// Section 2: pad counts in the XZ and YZ projections
// ============================================================================
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
inline PadProjectionSummary ComputePadProjectionCounts(const std::vector<ActRoot::Voxel>& voxels, double padPitch = 2.0)
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

// ============================================================================
// Section 3: gaps in a "row" coordinate per Z-slice, restricted to a window
// at the END of the track (largest s), where the Bragg peak deposits most of
// the charge.
// ============================================================================
//
// GENERIC: it doesn't hardcode X or Y. You pass in a lambda (getRow) that
// extracts whichever coordinate you want to use as the "row" grouped by Z
// (either pos.X() for the XZ projection or pos.Y() for the YZ projection).
// For each fixed z (a "row" of pads at a given drift time) *within that
// window*, take the sorted, distinct row-coordinate values that are occupied.
// Any place where two consecutive occupied values are not adjacent (diff > 1
// pad unit) counts as ONE hole, no matter how many pads are missing in
// between.
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

inline ProjGapResult ComputeGapsProjectionInRange(const std::vector<ActRoot::Voxel>& voxels,
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

// --- choose, per event, whichever projection (XZ or YZ) is "wide" ---------
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

inline ChosenProjectionResult ComputeGapsDenserProjection(const std::vector<ActRoot::Voxel>& voxels,
                                                          double fracFromEnd = 0.5, int minDeltaZ = 0,
                                                          int maxHoleSize = 2)
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

#endif // GapAnalysisHelpers_h