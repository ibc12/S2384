#ifndef Check_dp_rpNearEnd_cxx
#define Check_dp_rpNearEnd_cxx

#include "ActDataManager.h"
#include "ActInputParser.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActSilData.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TString.h"
#include "TStyle.h"

#include "Math/Vector3D.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "../Multiplicity4/Utils.h"

namespace CheckDpRp
{
// ---------------- Configuration ----------------
constexpr double kPadToMm = 2.0;                        // 240 mm <-> 120 pad units
constexpr double kRpMinPad = 90.;                       // RP.X() > 120 pad units (= 240 mm)
constexpr double kRpMinMm = kRpMinPad * kPadToMm;       // RP.X() > 240 mm
constexpr double kMinMmPlot = kRpMinMm - 10;            // RP.X() < 120 pad units (= 240 mm)
constexpr double kMaxMmPlot = 260.;                     // max X in mm for the plots
constexpr double kThetaMinDeg = 10.;                    // min angle between light and beam
constexpr double kThetaMaxDegPlot = 185.;               // max angle between light and beam
constexpr double kThetaMinDegPlot = kThetaMinDeg - 10.; // max angle between light and beam
constexpr double kValidationZone = 4.;                  // pad units from the TPC edges (SET YOUR VALUE)
constexpr double kQaveMin = 400.;                       // Qave cut for L1
// If true, a GATCONF 1/2 event is "lateral" only if it does NOT stop inside the TPC
constexpr bool kLateralRequireNotStopped = true;
// -----------------------------------------------

constexpr double kRadToDeg = 180. / M_PI;

template <typename T>
ROOT::Math::XYZVector UnitDir(const T& v)
{
    ROOT::Math::XYZVector d {v.X(), v.Y(), v.Z()};
    if(d.Mag2() <= 0.)
        return d;
    return d.Unit();
}

// Beam = cluster whose direction is closest to (1,0,0)
int GetBeamIdx(ActRoot::TPCData& tpc)
{
    int idx {-1};
    double best {-1.};
    for(int i = 0; i < (int)tpc.fClusters.size(); ++i)
    {
        auto d {UnitDir(tpc.fClusters[i].GetLine().GetDirection())};
        double c {std::abs(d.X())};
        if(c > best)
        {
            best = c;
            idx = i;
        }
    }
    return idx;
}

// Direction of a cluster oriented from the RP towards its farthest voxel
ROOT::Math::XYZVector GetOrientedDir(ActRoot::TPCData& tpc, int idx)
{
    auto& cluster {tpc.fClusters.at(idx)};
    auto dir {UnitDir(cluster.GetLine().GetDirection())};
    auto rp {tpc.fRPs.front()};
    double maxDist {-1.};
    ROOT::Math::XYZVector far {};
    for(auto& v : cluster.GetRefToVoxels())
    {
        auto p {v.GetPosition()};
        ROOT::Math::XYZVector d {p.X() - rp.X(), p.Y() - rp.Y(), p.Z() - rp.Z()};
        if(d.Mag2() > maxDist)
        {
            maxDist = d.Mag2();
            far = d;
        }
    }
    if(dir.Dot(far) < 0.)
        dir = -dir;
    return dir;
}

double ToDeg(double cosine)
{
    return std::acos(std::clamp(cosine, -1., 1.)) * kRadToDeg;
}

// Geometric check only: the last voxel of the light cluster (sorted along its direction)
// must be inside the TPC limits (minus the validation zone) in X and Y.
// The direction is oriented from the RP so that voxels.back() is the end of the track.
bool StopsInsideTPC(ActRoot::TPCData& tpc, int lightIdx)
{
    if(lightIdx < 0)
        return false;
    auto cluster = tpc.fClusters[lightIdx]; // local copy, does not modify the original
    auto dir = cluster.GetLine().GetDirection();

    if(!tpc.fRPs.empty())
    {
        auto rp {tpc.fRPs.front()};
        double maxDist {-1.};
        double dot {0.};
        for(auto& v : cluster.GetRefToVoxels())
        {
            auto p {v.GetPosition()};
            double dx {p.X() - rp.X()}, dy {p.Y() - rp.Y()}, dz {p.Z() - rp.Z()};
            double d2 {dx * dx + dy * dy + dz * dz};
            if(d2 > maxDist)
            {
                maxDist = d2;
                dot = dir.X() * dx + dir.Y() * dy + dir.Z() * dz;
            }
        }
        if(dot < 0.)
            dir = -dir;
    }

    cluster.SortAlongDir(dir);
    auto& voxels = cluster.GetRefToVoxels();
    if(voxels.empty())
        return false;
    auto stopPoint = voxels.back().GetPosition();
    bool insideX = (stopPoint.X() > kValidationZone && stopPoint.X() < (128 - kValidationZone));
    bool insideY = (stopPoint.Y() > kValidationZone && stopPoint.Y() < (128 - kValidationZone));
    return insideX && insideY;
}

// Silicon layer associated to each GATCONF (empty if none)
std::string LayerOfGatconf(int g)
{
    if(g == 1)
        return "l0";
    if(g == 2)
        return "r0";
    return "";
}
} // namespace CheckDpRp

void Check_dp_rpNearEnd_Mult2()
{
    using namespace CheckDpRp;

    std::string beam {"7Li"};

    std::string dataconf {};
    if(beam == "11Li")
        dataconf = "./../../configs/data_11Li.conf";
    else if(beam == "7Li")
        dataconf = "./../../configs/data_7Li.conf";
    else
        throw std::runtime_error("Beam cannot differ from 11Li or 7Li");

    // Read data: filter (clusters + RP), silicon + modular, merger and raw TPC data
    ActRoot::DataManager dataman {dataconf, ActRoot::ModeType::EFilter};
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    auto chain3 {dataman.GetChain(ActRoot::ModeType::EMerge)};
    auto chain4 {dataman.GetChain(ActRoot::ModeType::EReadTPC)};
    chain->AddFriend(chain2.get());
    chain->AddFriend(chain3.get());
    chain->AddFriend(chain4.get(), "GETTree");

    // Get drift parameter
    ActRoot::InputParser parser {};
    parser.ReadFile("../../configs/detector.conf");
    auto mergerBlock = parser.GetBlock("Merger");
    auto driftFactor = mergerBlock->GetDouble("DriftFactor"); // in pad units

    // Read per-run zDrift cuts: columns are run, meanStart, sigmaStart, meanEnd, sigmaEnd
    struct ZDriftCut
    {
        double zMin; // meanStart + nSigma * sigmaStart
        double zMax; // meanEnd   - nSigma * sigmaEnd
    };
    const double nSigmaZDrift {3.0};
    std::map<int, ZDriftCut> zDriftCuts;
    {
        std::ifstream finZ("../L1PID/Outputs/zDrift_cut_perRun.dat");
        if(!finZ.is_open())
            throw std::runtime_error("Could not open zDrift_cut_perRun.dat");
        std::string line;
        while(std::getline(finZ, line))
        {
            if(line.empty() || line[0] == '#')
                continue;
            std::istringstream iss(line);
            int run;
            double meanS, sigmaS, meanE, sigmaE;
            if(iss >> run >> meanS >> sigmaS >> meanE >> sigmaE)
                zDriftCuts[run] = {meanS + nSigmaZDrift * sigmaS, meanE - nSigmaZDrift * sigmaE};
        }
        finZ.close();
    }

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {*chain};

    auto dfSel =
        df.Filter("fClusters.size() == 2")
            // Only GATCONF 1 (l0), 2 (r0) and 8 (L1, stopped inside)
            .Filter(
                [](ActRoot::ModularData& m)
                {
                    int g {static_cast<int>(std::lround(m.Get("GATCONF")))};
                    return g == 1 || g == 2 || g == 8;
                },
                {"ModularData"})
            .Define("GATCONF", [](ActRoot::ModularData& m) { return static_cast<int>(std::lround(m.Get("GATCONF"))); },
                    {"ModularData"})
            // Remove events with high charge deposit around the RP
            // (most of the times it is masked by rp or is in beam cluster)
            .Filter(
                [&](ActRoot::TPCData& f, ActRoot::TPCData& tpc)
                {
                    // Check first if the RP is defined, otherwise return false
                    if(f.fRPs.empty())
                        return false;
                    auto rp {f.fRPs.front()};
                    auto rp_x {rp.X()};
                    auto rp_y {rp.Y()};
                    // Run for all clusters
                    int counter {};
                    for(auto& cluster : tpc.fClusters)
                    {
                        auto voxels {cluster.GetRefToVoxels()};
                        for(auto& v : voxels)
                        {
                            if(v.GetPosition().Y() > rp_y - 5 && v.GetPosition().Y() < rp_y + 5 &&
                               v.GetPosition().X() > rp_x - 10 && v.GetPosition().X() < rp_x + 10)
                                if(v.GetCharge() > 3000.)
                                    counter++;
                        }
                    }
                    if(counter > 4)
                        return false;
                    return true;
                },
                {"TPCData", "GETTree.TPCData"})
            // RP defined and RP.X() > 120 pad units
            .Filter([](ActRoot::TPCData& tpc) { return !tpc.fRPs.empty(); }, {"TPCData"})
            .Define("RPx", [](ActRoot::TPCData& tpc) { return (double)tpc.fRPs.front().X(); }, {"TPCData"})
            .Filter([](double rpx) { return rpx > kRpMinPad; }, {"RPx"})
            .Define("RPxmm", [](double rpx) { return rpx * kPadToMm; }, {"RPx"})
            // Beam-like = closest to (1,0,0); light = the other one
            .Define("BeamIdx", GetBeamIdx, {"TPCData"})
            .Filter([](int b) { return b >= 0; }, {"BeamIdx"})
            .Define("LightIdx", [](int b) { return 1 - b; }, {"BeamIdx"})
            // Angle between light and beam directions
            .Define("ThetaBeam",
                    [](ActRoot::TPCData& tpc, int b, int l)
                    {
                        auto db {GetOrientedDir(tpc, b)};
                        if(db.X() < 0.) // beam goes towards +X
                            db = -db;
                        auto dl {GetOrientedDir(tpc, l)};
                        return ToDeg(db.Dot(dl));
                    },
                    {"TPCData", "BeamIdx", "LightIdx"})
            .Filter([](double th) { return th > kThetaMinDeg; }, {"ThetaBeam"})
            // Angle of the light particle with respect to the X axis
            .Define("ThetaX",
                    [](ActRoot::TPCData& tpc, int l)
                    {
                        auto dl {GetOrientedDir(tpc, l)};
                        return ToDeg(dl.X());
                    },
                    {"TPCData", "LightIdx"})
            // zDrift of the light particle (RP forced to be at 110 mm)
            .Define("zDrift",
                    [&driftFactor](ActRoot::TPCData& tpc, int idx)
                    {
                        auto& cluster = tpc.fClusters.at(idx);
                        auto& voxels = cluster.GetRefToVoxels();

                        if(tpc.fRPs.empty())
                            return -111.0f;
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

                        float zDrift = 110.0 + deltaZ; // Force rp to be at 110 mm

                        return zDrift;
                    },
                    {"TPCData", "LightIdx"})
            // Average charge per unit length of the light particle
            .Define("Qave",
                    [](ActRoot::TPCData& tpc, int idx)
                    {
                        auto cluster = tpc.fClusters[idx];
                        cluster.SortAlongDir(cluster.GetLine().GetDirection());
                        auto voxels = cluster.GetRefToVoxels();
                        if(voxels.empty())
                            return -1.f;
                        float Q = 0.f;
                        for(const auto& v : voxels)
                            Q += v.GetCharge();
                        ROOT::Math::XYZPointF firstPos = voxels.front().GetPosition();
                        ROOT::Math::XYZPointF lastPos = voxels.back().GetPosition();
                        Utils::ScalePoint(firstPos);
                        Utils::ScalePoint(lastPos);
                        float length = (lastPos - firstPos).R();
                        if(length <= 0.0)
                            return -1.f;
                        return Q / length;
                    },
                    {"TPCData", "LightIdx"})
            // Geometric check: does the light particle stop inside the TPC limits?
            .Define("StopsInsideTPC", [](ActRoot::TPCData& tpc, int l) { return StopsInsideTPC(tpc, l); },
                    {"TPCData", "LightIdx"})
            // L1 condition: stops inside the TPC + zDrift and Qave inside the cuts per run
            .Define("StoppedInACTAR",
                    [&](ActRoot::MergerData& m, bool stopsInside, float zDrift, float Qave)
                    {
                        if(!stopsInside)
                            return false;

                        // Now check if the zDrift and Qave enters the cuts per run
                        bool isInsideZDriftCut = false;
                        bool isInsideQaveCut = false;
                        auto it = zDriftCuts.find(m.fRun);
                        if(it == zDriftCuts.end())
                            isInsideZDriftCut = false; // run not in file -> reject
                        else
                            isInsideZDriftCut = zDrift > it->second.zMin && zDrift < it->second.zMax;

                        isInsideQaveCut = Qave > kQaveMin;

                        return isInsideZDriftCut && isInsideQaveCut;
                    },
                    {"MergerData", "StopsInsideTPC", "zDrift", "Qave"})
            // Lateral: GATCONF 1 (l0) or 2 (r0) + hit in the corresponding silicon layer
            // (+ not stopping inside the TPC, if requested; geometric check only)
            .Define("isLateral",
                    [](ActRoot::SilData& sil, int g, bool stopsInside)
                    {
                        auto layer {LayerOfGatconf(g)};
                        if(layer.empty())
                            return false;
                        auto it {sil.fSiN.find(layer)};
                        bool hit {it != sil.fSiN.end() && !it->second.empty()};
                        if(!hit)
                            return false;
                        if(kLateralRequireNotStopped && stopsInside)
                            return false;
                        return true;
                    },
                    {"SilData", "GATCONF", "StopsInsideTPC"})
            // L1: GATCONF 8 + StoppedInACTAR
            .Define("isStopped", [](int g, bool stopped) { return g == 8 && stopped; }, {"GATCONF", "StoppedInACTAR"});

    auto dfLat {dfSel.Filter("isLateral")};
    auto dfStop {dfSel.Filter("isStopped")};

    // Counts
    auto nSel {dfSel.Count()};
    auto nG1 {dfSel.Filter("GATCONF == 1").Count()};
    auto nG2 {dfSel.Filter("GATCONF == 2").Count()};
    auto nG8 {dfSel.Filter("GATCONF == 8").Count()};
    auto nLat {dfLat.Count()};
    auto nStop {dfStop.Count()};
    auto nLatStopsInside {dfSel.Filter("(GATCONF == 1 || GATCONF == 2) && StopsInsideTPC").Count()};
    auto nG8NotStopped {dfSel.Filter("GATCONF == 8 && !StoppedInACTAR").Count()};

    // Histograms
    auto hTheta2D {dfSel.Histo2D({"hTheta2D", "#theta_{x} vs RP;RP.X [mm];#theta_{x} [#circ]", 68, kMinMmPlot,
                                  kMaxMmPlot, 100, kThetaMinDegPlot, kThetaMaxDegPlot},
                                 "RPxmm", "ThetaX")};
    auto hTheta {dfSel.Histo1D(
        {"hTheta", "#theta_{x};#theta_{x} [#circ];Counts", 100, kThetaMinDegPlot, kThetaMaxDegPlot}, "ThetaX")};
    auto hThetaLat {dfLat.Histo1D({"hThetaLat", "", 100, kThetaMinDegPlot, kThetaMaxDegPlot}, "ThetaX")};
    auto hThetaStop {dfStop.Histo1D({"hThetaStop", "", 100, kThetaMinDegPlot, kThetaMaxDegPlot}, "ThetaX")};
    auto hRP {dfSel.Histo1D({"hRP", "RP;RP.X [mm];Counts", 68, kMinMmPlot, kMaxMmPlot}, "RPxmm")};
    auto hRPLat {dfLat.Histo1D({"hRPLat", "", 68, kMinMmPlot, kMaxMmPlot}, "RPxmm")};
    auto hRPStop {dfStop.Histo1D({"hRPStop", "", 68, kMinMmPlot, kMaxMmPlot}, "RPxmm")};
    auto h2Lat {dfLat.Histo2D({"h2Lat", "Lateral;RP.X [mm];#theta_{x} [#circ]", 68, kMinMmPlot, kMaxMmPlot, 100,
                               kThetaMinDegPlot, kThetaMaxDegPlot},
                              "RPxmm", "ThetaX")};
    auto h2Stop {dfStop.Histo2D({"h2Stop", "Stopped;RP.X [mm];#theta_{x} [#circ]", 68, kMinMmPlot, kMaxMmPlot, 100,
                                 kThetaMinDegPlot, kThetaMaxDegPlot},
                                "RPxmm", "ThetaX")};

    // Print counts (this triggers the event loop for everything booked above)
    std::cout << "---- Check_dp_rpNearEnd (" << beam << ") ----" << '\n';
    std::cout << "Selected events (M2, GATCONF 1/2/8, charge filter, RP.X > " << kRpMinPad << " pad, theta > "
              << kThetaMinDeg << " deg): " << nSel.GetValue() << '\n';
    std::cout << "  GATCONF 1 (l0): " << nG1.GetValue() << '\n';
    std::cout << "  GATCONF 2 (r0): " << nG2.GetValue() << '\n';
    std::cout << "  GATCONF 8 (L1): " << nG8.GetValue() << '\n';
    std::cout << "-> Lateral (hit in l0/r0" << (kLateralRequireNotStopped ? ", not stopping inside" : "")
              << "): " << nLat.GetValue() << '\n';
    std::cout << "-> L1 stopped inside (GATCONF 8 + StoppedInACTAR): " << nStop.GetValue() << '\n';
    std::cout << "(GATCONF 1/2 that stop inside the TPC: " << nLatStopsInside.GetValue()
              << " | GATCONF 8 that do NOT pass StoppedInACTAR: " << nG8NotStopped.GetValue() << ")" << '\n';

    // Styles
    hTheta->SetLineColor(kBlack);
    hThetaLat->SetLineColor(kRed);
    hThetaStop->SetLineColor(kBlue);
    hRP->SetLineColor(kBlack);
    hRPLat->SetLineColor(kRed);
    hRPStop->SetLineColor(kBlue);
    h2Lat->SetMarkerColor(kRed);
    h2Stop->SetMarkerColor(kBlue);
    h2Lat->SetMarkerStyle(20);
    h2Stop->SetMarkerStyle(20);
    h2Lat->SetMarkerSize(0.6);
    h2Stop->SetMarkerSize(0.6);

    // Plots
    auto* c = new TCanvas("cCheckDp", "dp RP near end", 1200, 900);
    c->Divide(2, 2);

    c->cd(1);
    hTheta2D->DrawClone("colz");

    c->cd(2);
    auto* leg = new TLegend(0.6, 0.7, 0.88, 0.88);
    hTheta->DrawClone("hist");
    hThetaLat->DrawClone("hist same");
    hThetaStop->DrawClone("hist same");
    leg->AddEntry(hTheta->Clone(), "All", "l");
    leg->AddEntry(hThetaLat->Clone(), "Lateral (l0/r0)", "l");
    leg->AddEntry(hThetaStop->Clone(), "L1 stopped inside", "l");
    leg->Draw();

    c->cd(3);
    hRP->DrawClone("hist");
    hRPLat->DrawClone("hist same");
    hRPStop->DrawClone("hist same");

    c->cd(4);
    auto* frame = (TH2D*)h2Lat->Clone("frame");
    frame->Reset();
    frame->SetTitle("Lateral (red) vs L1 stopped (blue);RP.X [mm];#theta_{x} [#circ]");
    frame->SetStats(false);
    frame->DrawClone("axis");
    h2Lat->DrawClone("scat same");
    h2Stop->DrawClone("scat same");

    // Save events
    // std::ofstream out(TString::Format("./Outputs/events_mult2_rpNearEnd_%s_dp_lat.dat", beam.c_str()));
    // dfLat.Foreach([&](ActRoot::MergerData& m) { m.Stream(out); }, {"MergerData"});
    // out.close();
    // std::ofstream out1(TString::Format("./Outputs/events_mult2_rpNearEnd_%s_dp_L1.dat", beam.c_str()));
    // dfStop.Foreach([&](ActRoot::MergerData& m) { m.Stream(out1); }, {"MergerData"});
    // out1.close();
}
#endif