#ifndef CHECKRAWTL_TL_DEFINITIONS_H
#define CHECKRAWTL_TL_DEFINITIONS_H

#include "ActCluster.h"
#include "ActDataManager.h"
#include "ActInputParser.h"
#include "ActLine.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TH2.h"
#include "TStyle.h"

#include "Math/Point3D.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace CheckTL
{
using XYZPointF = ROOT::Math::XYZPointF;

// Same as ScalePoint(p, side, drift, true) in MergerDetector: voxel center offset + scale
XYZPointF ScaleWithOffset(const XYZPointF& p, float padSide, float driftFactor)
{
    return {(p.X() + 0.5f) * padSide, (p.Y() + 0.5f) * padSide, (p.Z() + 0.5f) * driftFactor};
}

// Track length of a cluster projected onto its fitted line.
//   beginAtRP = false -> begin = first voxel of cluster   (original RawTL definition)
//   beginAtRP = true  -> begin = reaction point           (original TL definition)
//   scale     = false -> pad units (offset correction only)
//   scale     = true  -> mm (x, y * padSide; z * driftFactor)
// rpScaled is MergerData::fRP, which is already in mm
double ComputeTL(ActRoot::Cluster cluster, // by value: we sort a copy
                 const XYZPointF& rpScaled, bool beginAtRP, bool scale, float padSide, float driftFactor)
{
    auto& voxels {cluster.GetRefToVoxels()};
    if(voxels.empty())
        return -1.;
    // Sorting is done in LightOrHeavy in the Merger, but the TPCData tree may be stored before it
    cluster.SortAlongDir();

    XYZPointF begin {};
    if(beginAtRP)
    {
        if(scale)
            begin = rpScaled;
        else // undo the scaling to get back to pad units (offset already included)
            begin = {rpScaled.X() / padSide, rpScaled.Y() / padSide, rpScaled.Z() / driftFactor};
    }
    else
    {
        begin = voxels.front().GetPosition();
        begin = scale ? ScaleWithOffset(begin, padSide, driftFactor) : ScaleWithOffset(begin, 1, 1);
    }
    auto end {voxels.back().GetPosition()};
    end = scale ? ScaleWithOffset(end, padSide, driftFactor) : ScaleWithOffset(end, 1, 1);

    auto line {cluster.GetLine()}; // copy
    if(scale)
        line.Scale(padSide, driftFactor);

    auto projBegin {line.ProjectionPointOnLine(begin)};
    auto projEnd {line.ProjectionPointOnLine(end)};
    return (projBegin - projEnd).R();
}
} // namespace CheckTL

void CheckRawTL_TL_Definitions(const std::string& beam = "7Li", bool isLight = true)
{
    auto inFile {TString::Format("../../PostAnalysis/Outputs/tree_preprocess_F_%s.root", beam.c_str())};

    // RDataFrame
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame dfin {"PreProcessed_Tree", inFile};
    auto df = dfin.Filter([=](ActRoot::ModularData& m, ActRoot::MergerData& mer) { return m.Get("GATCONF") == 8 && (mer.fLightIdx != -1); },
                      {"ModularData", "MergerData"});

    // Scale factors
    const float padSide {2.f};
    ActRoot::InputParser parser {};
    parser.ReadFile("../../configs/detector.conf");
    // TODO: replace by however you end up reading the drift factor
    const float driftFactor {static_cast<float>(parser.GetBlock("Merger")->GetDouble("DriftFactor"))};


    // Compute the 4 track lengths at once: {RawTL_old, TL_new, RawTL_new, TL_old}
    //   [0] first voxel -> last voxel, unscaled  (original RawTL)
    //   [1] first voxel -> last voxel, scaled    (new TL)
    //   [2] RP          -> last voxel, unscaled  (new RawTL)
    //   [3] RP          -> last voxel, scaled    (original TL)
    auto dfTL =
        df.Define("TLs",
                  [=](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                  {
                      std::vector<double> out(4, -1.);
                      int idx {isLight ? m.fLightIdx : m.fHeavyIdx};
                      if(idx < 0 || static_cast<std::size_t>(idx) >= tpc.fClusters.size())
                          return out;
                      auto& cl {tpc.fClusters.at(idx)};
                      out[0] = CheckTL::ComputeTL(cl, m.fRP, false, false, padSide, driftFactor);
                      out[1] = CheckTL::ComputeTL(cl, m.fRP, false, true, padSide, driftFactor);
                      out[2] = CheckTL::ComputeTL(cl, m.fRP, true, false, padSide, driftFactor);
                      out[3] = CheckTL::ComputeTL(cl, m.fRP, true, true, padSide, driftFactor);
                      return out;
                  },
                  {"MergerData", "TPCData"})
            .Filter([](const std::vector<double>& v) { return v[0] > 0 && v[1] > 0 && v[2] > 0 && v[3] > 0; }, {"TLs"})
            .Define("Qtotal", [=](ActRoot::MergerData& m) -> double
                    { return isLight ? m.fLight.fQtotal : m.fHeavy.fQtotal; }, {"MergerData"})
            .Define("StoredRawTL", [=](ActRoot::MergerData& m) -> double
                    { return isLight ? m.fLight.fRawTL : m.fHeavy.fRawTL; }, {"MergerData"})
            .Define("StoredTL", [=](ActRoot::MergerData& m) -> double { return isLight ? m.fLight.fTL : m.fHeavy.fTL; },
                    {"MergerData"})
            .Define("RawTL_old", [](const std::vector<double>& v) { return v[0]; }, {"TLs"})
            .Define("TL_new", [](const std::vector<double>& v) { return v[1]; }, {"TLs"})
            .Define("RawTL_new", [](const std::vector<double>& v) { return v[2]; }, {"TLs"})
            .Define("TL_old", [](const std::vector<double>& v) { return v[3]; }, {"TLs"});

    // Sanity check: do my "old" definitions reproduce what is stored in MergerData?
    auto nTotal = dfTL.Count();
    auto nBadRaw = dfTL.Filter([](double mine, double stored) { return std::abs(mine - stored) > 1e-2; },
                               {"RawTL_old", "StoredRawTL"})
                       .Count();
    auto nBadTL =
        dfTL.Filter([](double mine, double stored) { return std::abs(mine - stored) > 1e-2; }, {"TL_old", "StoredTL"})
            .Count();

    // Histograms (edit ranges to taste)
    const int nBinsTL {300}, nBinsQ {300};
    const double tlMaxPad {150.}, tlMaxMM {300.}, qMax {60000.};
    auto hRawOld =
        dfTL.Histo2D({"hRawOld", "ORIGINAL RawTL: first voxel #rightarrow last, unscaled;RawTL [pad units];Qtotal",
                      200, 0, 120, 2000, 0, 3e5},
                     "RawTL_old", "Qtotal");
    auto hTLNew = dfTL.Histo2D(
        {"hTLNew", "NEW TL: first voxel #rightarrow last, scaled;TL [mm];Qtotal", 300, 0, 200, 2000, 0, 3e5},
        "TL_new", "Qtotal");
    auto hRawNew = dfTL.Histo2D({"hRawNew", "NEW RawTL: RP #rightarrow last, unscaled;RawTL [pad units];Qtotal",
                                 200, 0, 120, 2000, 0, 3e5},
                                "RawTL_new", "Qtotal");
    auto hTLOld = dfTL.Histo2D(
        {"hTLOld", "ORIGINAL TL: RP #rightarrow last, scaled;TL [mm];Qtotal", 300, 0, 200, 2000, 0, 3e5},
        "TL_old", "Qtotal");

    // Triggers the event loop (once for everything)
    std::cout << "Events used: " << nTotal.GetValue() << '\n';
    std::cout << "Mismatch vs stored fRawTL: " << nBadRaw.GetValue() << '\n';
    std::cout << "Mismatch vs stored fTL   : " << nBadTL.GetValue() << '\n';

    gStyle->SetOptStat(0);
    auto* c {new TCanvas {"cTL", "Qtotal vs track length definitions", 1400, 1000}};
    c->Divide(2, 2);
    c->cd(1);
    gPad->SetLogz();
    hRawOld->DrawClone("colz");
    c->cd(2);
    gPad->SetLogz();
    hTLNew->DrawClone("colz");
    c->cd(3);
    gPad->SetLogz();
    hRawNew->DrawClone("colz");
    c->cd(4);
    gPad->SetLogz();
    hTLOld->DrawClone("colz");
}

#endif