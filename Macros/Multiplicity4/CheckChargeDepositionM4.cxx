#include "ActCluster.h"
#include "ActMergerData.h"
#include "ActTPCData.h"
#include "ActVoxel.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"

#include <fstream>
#include <iostream>

#include "../../PrettyStyle.C"


void ScalePoint(ROOT::Math::XYZPointF& point, float xy, float z, bool addOffset = false)
{
    if(addOffset) // when converting a bin point to physical units which wasnt already corrected
        point += ROOT::Math::XYZVector {0.5, 0.5, 0.5};
    point.SetX(point.X() * xy);
    point.SetY(point.Y() * xy);
    point.SetZ(point.Z() * z);
}

void CheckChargeDepositionM4()
{
    PrettyStyle(true);
    // ROOT::DisableImplicitMT();
    ROOT::RDataFrame df {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_7Li_d_p_filtered.root"};

    // Start with only silicon events and ensure M=4
    auto dfFilter =
        df.Filter([](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                  { return (m.fLight.IsFilled() == true && tpc.fClusters.size() == 4); }, {"MergerData", "TPCData"});

    //  Define the charge / dist for each particle except light particle detected in silicon

    auto dfIndex = dfFilter
                       .Define("MaxThetaIdx",
                               [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                               {
                                   int lightIdx = m.fLightIdx;
                                   int beamIdx = m.fBeamIdx;
                                   auto clusters = tpc.fClusters;
                                   int idx = -1;
                                   double maxAngle = -1.0;
                                   for(int i = 0; i < clusters.size(); ++i)
                                   {
                                       if(i == lightIdx || i == beamIdx)
                                           continue;
                                       auto beamDir = clusters[beamIdx].GetLine().GetDirection().Unit();
                                       auto currentDir = clusters[i].GetLine().GetDirection().Unit();
                                       double angle = TMath::ACos(beamDir.Dot(currentDir)) * TMath::RadToDeg();
                                       if(angle > maxAngle)
                                       {
                                           maxAngle = angle;
                                           idx = i;
                                       }
                                   }
                                   return idx;
                               },
                               {"MergerData", "TPCData"})
                       .Define("MinThetaIdx",
                               [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                               {
                                   int lightIdx = m.fLightIdx;
                                   int beamIdx = m.fBeamIdx;
                                   auto clusters = tpc.fClusters;
                                   int idx = -1;
                                   double minAngle = 1e9;
                                   for(int i = 0; i < clusters.size(); ++i)
                                   {
                                       if(i == lightIdx || i == beamIdx)
                                           continue;
                                       auto beamDir = clusters[beamIdx].GetLine().GetDirection().Unit();
                                       auto currentDir = clusters[i].GetLine().GetDirection().Unit();
                                       double angle = TMath::ACos(beamDir.Dot(currentDir)) * TMath::RadToDeg();
                                       if(angle < minAngle)
                                       {
                                           minAngle = angle;
                                           idx = i;
                                       }
                                   }
                                   return idx;
                               },
                               {"MergerData", "TPCData"})
                       .Define("MaxQLengthIdx",
                               [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                               {
                                   int lightIdx = m.fLightIdx;
                                   int beamIdx = m.fBeamIdx;
                                   auto clusters = tpc.fClusters;
                                   int idx = -1;
                                   double qlMax = -1.0;
                                   for(int i = 0; i < clusters.size(); ++i)
                                   {
                                       if(i == lightIdx || i == beamIdx)
                                           continue;
                                       auto& cl = clusters[i];
                                       cl.SortAlongDir(cl.GetLine().GetDirection());
                                       auto voxels = cl.GetRefToVoxels();
                                       if(voxels.size() < 2)
                                           continue;
                                       float totalCharge = 0;
                                       for(const auto& v : voxels)
                                           totalCharge += v.GetCharge();
                                       ROOT::Math::XYZPointF firstPos = voxels.front().GetPosition();
                                       ROOT::Math::XYZPointF lastPos = voxels.back().GetPosition();
                                       ScalePoint(firstPos, 2, 2.84032);
                                       ScalePoint(lastPos, 2, 2.84032);
                                       double length = (lastPos - firstPos).R();
                                       if(length <= 0.0)
                                           continue;
                                       double ql = totalCharge / length;
                                       if(ql > qlMax)
                                       {
                                           qlMax = ql;
                                           idx = i;
                                       }
                                   }
                                   return idx;
                               },
                               {"MergerData", "TPCData"})
                       .Define("MinQLengthIdx",
                               [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                               {
                                   int lightIdx = m.fLightIdx;
                                   int beamIdx = m.fBeamIdx;
                                   auto clusters = tpc.fClusters;
                                   int idx = -1;
                                   double qlMin = 1e9;
                                   for(int i = 0; i < clusters.size(); ++i)
                                   {
                                       if(i == lightIdx || i == beamIdx)
                                           continue;
                                       auto& cl = clusters[i];
                                       cl.SortAlongDir(cl.GetLine().GetDirection());
                                       auto voxels = cl.GetRefToVoxels();
                                       if(voxels.size() < 2)
                                           continue;
                                       float totalCharge = 0;
                                       for(const auto& v : voxels)
                                           totalCharge += v.GetCharge();
                                       ROOT::Math::XYZPointF firstPos = voxels.front().GetPosition();
                                       ROOT::Math::XYZPointF lastPos = voxels.back().GetPosition();
                                       ScalePoint(firstPos, 2, 2.84032);
                                       ScalePoint(lastPos, 2, 2.84032);
                                       double length = (lastPos - firstPos).R();
                                       if(length <= 0.0)
                                           continue;
                                       double ql = totalCharge / length;
                                       if(ql < qlMin)
                                       {
                                           qlMin = ql;
                                           idx = i;
                                       }
                                   }
                                   return idx;
                               },
                               {"MergerData", "TPCData"});

    auto dfCharge = dfIndex
                        .Define("ChargeBeam",
                                [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                                {
                                    int beamIdx = m.fBeamIdx;
                                    auto& cluster = tpc.fClusters[beamIdx];
                                    cluster.SortAlongDir(cluster.GetLine().GetDirection());
                                    auto voxels = cluster.GetRefToVoxels();
                                    float totalCharge = 0.f;
                                    for(const auto& v : voxels)
                                        totalCharge += v.GetCharge();
                                    ROOT::Math::XYZPointF firstPos = voxels.front().GetPosition();
                                    ROOT::Math::XYZPointF lastPos = voxels.back().GetPosition();
                                    ScalePoint(firstPos, 2, 2.84032);
                                    ScalePoint(lastPos, 2, 2.84032);
                                    float length = (lastPos - firstPos).R();
                                    if(length <= 0.0)
                                        return 0.0f;
                                    return totalCharge / length;
                                },
                                {"MergerData", "TPCData"})
                        .Define("ChargeLight1",
                                [](ActRoot::MergerData& m, ActRoot::TPCData& tpc, int idx)
                                {
                                    if(idx < 0)
                                        return -1.f;
                                    auto& cluster = tpc.fClusters[idx];
                                    cluster.SortAlongDir(cluster.GetLine().GetDirection());
                                    auto voxels = cluster.GetRefToVoxels();
                                    float totalCharge = 0.f;
                                    for(const auto& v : voxels)
                                        totalCharge += v.GetCharge();
                                    ROOT::Math::XYZPointF firstPos = voxels.front().GetPosition();
                                    ROOT::Math::XYZPointF lastPos = voxels.back().GetPosition();
                                    ScalePoint(firstPos, 2, 2.84032);
                                    ScalePoint(lastPos, 2, 2.84032);
                                    float length = (lastPos - firstPos).R();
                                    if(length <= 0.0)
                                        return -1.f;
                                    return totalCharge / length;
                                },
                                {"MergerData", "TPCData", "MaxThetaIdx"})
                        .Define("ChargeLight2",
                                [](ActRoot::MergerData& m, ActRoot::TPCData& tpc, int idx)
                                {
                                    if(idx < 0)
                                        return -1.f;
                                    auto& cluster = tpc.fClusters[idx];
                                    cluster.SortAlongDir(cluster.GetLine().GetDirection());
                                    auto voxels = cluster.GetRefToVoxels();
                                    float totalCharge = 0.f;
                                    for(const auto& v : voxels)
                                        totalCharge += v.GetCharge();
                                    ROOT::Math::XYZPointF firstPos = voxels.front().GetPosition();
                                    ROOT::Math::XYZPointF lastPos = voxels.back().GetPosition();
                                    ScalePoint(firstPos, 2, 2.84032);
                                    ScalePoint(lastPos, 2, 2.84032);
                                    float length = (lastPos - firstPos).R();
                                    if(length <= 0.0)
                                        return -1.f;
                                    return totalCharge / length;
                                },
                                {"MergerData", "TPCData", "MinThetaIdx"})
                        .Define("ThetaLight1",
                                [](ActRoot::MergerData& m, ActRoot::TPCData& tpc, int idx)
                                {
                                    if(idx < 0)
                                        return -1.0;
                                    auto& cluster = tpc.fClusters[idx];
                                    auto beamDir = tpc.fClusters[m.fBeamIdx].GetLine().GetDirection().Unit();
                                    double cosang = beamDir.Dot(cluster.GetLine().GetDirection().Unit());
                                    return TMath::ACos(cosang) * TMath::RadToDeg();
                                },
                                {"MergerData", "TPCData", "MaxThetaIdx"})
                        .Define("ThetaLight2",
                                [](ActRoot::MergerData& m, ActRoot::TPCData& tpc, int idx)
                                {
                                    if(idx < 0)
                                        return -1.0;
                                    auto& cluster = tpc.fClusters[idx];
                                    auto beamDir = tpc.fClusters[m.fBeamIdx].GetLine().GetDirection().Unit();
                                    double cosang = beamDir.Dot(cluster.GetLine().GetDirection().Unit());
                                    return TMath::ACos(cosang) * TMath::RadToDeg();
                                },
                                {"MergerData", "TPCData", "MinThetaIdx"})
                        .Define("QLengthMax",
                                [](ActRoot::MergerData& m, ActRoot::TPCData& tpc, int idx)
                                {
                                    if(idx < 0)
                                        return -1.f;
                                    auto& cluster = tpc.fClusters[idx];
                                    cluster.SortAlongDir(cluster.GetLine().GetDirection());
                                    auto voxels = cluster.GetRefToVoxels();
                                    float totalCharge = 0.f;
                                    for(const auto& v : voxels)
                                        totalCharge += v.GetCharge();
                                    ROOT::Math::XYZPointF firstPos = voxels.front().GetPosition();
                                    ROOT::Math::XYZPointF lastPos = voxels.back().GetPosition();
                                    ScalePoint(firstPos, 2, 2.84032);
                                    ScalePoint(lastPos, 2, 2.84032);
                                    float length = (lastPos - firstPos).R();
                                    if(length <= 0.0)
                                        return -1.f;
                                    return totalCharge / length;
                                },
                                {"MergerData", "TPCData", "MaxQLengthIdx"})
                        .Define("QLengthMin",
                                [](ActRoot::MergerData& m, ActRoot::TPCData& tpc, int idx)
                                {
                                    if(idx < 0)
                                        return -1.f;
                                    auto& cluster = tpc.fClusters[idx];
                                    cluster.SortAlongDir(cluster.GetLine().GetDirection());
                                    auto voxels = cluster.GetRefToVoxels();
                                    float totalCharge = 0.f;
                                    for(const auto& v : voxels)
                                        totalCharge += v.GetCharge();
                                    ROOT::Math::XYZPointF firstPos = voxels.front().GetPosition();
                                    ROOT::Math::XYZPointF lastPos = voxels.back().GetPosition();
                                    ScalePoint(firstPos, 2, 2.84032);
                                    ScalePoint(lastPos, 2, 2.84032);
                                    float length = (lastPos - firstPos).R();
                                    if(length <= 0.0)
                                        return -1.f;
                                    return totalCharge / length;
                                },
                                {"MergerData", "TPCData", "MinQLengthIdx"})
                        .Define("ThetaQLengthMax",
                                [](ActRoot::MergerData& m, ActRoot::TPCData& tpc, int idx)
                                {
                                    if(idx < 0)
                                        return -1.0;
                                    auto& cluster = tpc.fClusters[idx];
                                    auto beamDir = tpc.fClusters[m.fBeamIdx].GetLine().GetDirection().Unit();
                                    double cosang = beamDir.Dot(cluster.GetLine().GetDirection().Unit());
                                    return TMath::ACos(cosang) * TMath::RadToDeg();
                                },
                                {"MergerData", "TPCData", "MaxQLengthIdx"})
                        .Define("ThetaQLengthMin",
                                [](ActRoot::MergerData& m, ActRoot::TPCData& tpc, int idx)
                                {
                                    if(idx < 0)
                                        return -1.0;
                                    auto& cluster = tpc.fClusters[idx];
                                    auto beamDir = tpc.fClusters[m.fBeamIdx].GetLine().GetDirection().Unit();
                                    double cosang = beamDir.Dot(cluster.GetLine().GetDirection().Unit());
                                    return TMath::ACos(cosang) * TMath::RadToDeg();
                                },
                                {"MergerData", "TPCData", "MinQLengthIdx"})
                        .Define("ChargeBeamAveQfrac",
                                [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                                {
                                    int beamIdx = m.fBeamIdx;
                                    auto& cluster = tpc.fClusters[beamIdx];
                                    auto voxels = cluster.GetRefToVoxels();
                                    if(voxels.size() < 5)
                                        return -100.f;
                                    auto dir = cluster.GetLine().GetDirection().Unit();
                                    ROOT::Math::XYZPointF r0 = voxels.front().GetPosition();
                                    std::vector<std::pair<float, float>> s_q;
                                    s_q.reserve(voxels.size());
                                    for(const auto& v : voxels)
                                        s_q.emplace_back((v.GetPosition() - r0).Dot(dir), v.GetCharge());
                                    float sMin = 1e9f, sMax = -1e9f;
                                    for(auto [s, q] : s_q)
                                    {
                                        sMin = std::min(sMin, s);
                                        sMax = std::max(sMax, s);
                                    }
                                    float L = sMax - sMin;
                                    if(L <= 0.f)
                                        return -100.f;
                                    float Qtotal = 0.f, Qtail = 0.f;
                                    float sCut = sMin + 0.7f * L;
                                    for(auto [s, q] : s_q)
                                    {
                                        Qtotal += q;
                                        if(s >= sCut)
                                            Qtail += q;
                                    }
                                    if(Qtotal <= 0.f)
                                        return -100.f;
                                    return Qtotal / Qtail;
                                },
                                {"MergerData", "TPCData"})
                        .Define("ChargeLight1AveQfrac",
                                [](ActRoot::MergerData& m, ActRoot::TPCData& tpc, int idx)
                                {
                                    if(idx < 0)
                                        return -1.f;
                                    auto& cluster = tpc.fClusters[idx];
                                    auto voxels = cluster.GetRefToVoxels();
                                    if(voxels.size() < 5)
                                        return -1.f;
                                    auto dir = cluster.GetLine().GetDirection().Unit();
                                    ROOT::Math::XYZPointF r0 = voxels.front().GetPosition();
                                    std::vector<std::pair<float, float>> s_q;
                                    s_q.reserve(voxels.size());
                                    for(const auto& v : voxels)
                                        s_q.emplace_back((v.GetPosition() - r0).Dot(dir), v.GetCharge());
                                    float sMin = 1e9f, sMax = -1e9f;
                                    for(auto [s, q] : s_q)
                                    {
                                        sMin = std::min(sMin, s);
                                        sMax = std::max(sMax, s);
                                    }
                                    float L = sMax - sMin;
                                    if(L <= 0.f)
                                        return -1.f;
                                    float Qtotal = 0.f, Qtail = 0.f;
                                    float sCut = sMin + 0.8f * L;
                                    for(auto [s, q] : s_q)
                                    {
                                        Qtotal += q;
                                        if(s >= sCut)
                                            Qtail += q;
                                    }
                                    if(Qtotal <= 0.f)
                                        return -1.f;
                                    return Qtotal / Qtail;
                                },
                                {"MergerData", "TPCData", "MaxThetaIdx"})
                        .Define("ChargeLight2AveQfrac",
                                [](ActRoot::MergerData& m, ActRoot::TPCData& tpc, int idx)
                                {
                                    if(idx < 0)
                                        return -1.f;
                                    auto& cluster = tpc.fClusters[idx];
                                    auto voxels = cluster.GetRefToVoxels();
                                    if(voxels.size() < 5)
                                        return -1.f;
                                    auto dir = cluster.GetLine().GetDirection().Unit();
                                    ROOT::Math::XYZPointF r0 = voxels.front().GetPosition();
                                    std::vector<std::pair<float, float>> s_q;
                                    s_q.reserve(voxels.size());
                                    for(const auto& v : voxels)
                                        s_q.emplace_back((v.GetPosition() - r0).Dot(dir), v.GetCharge());
                                    float sMin = 1e9f, sMax = -1e9f;
                                    for(auto [s, q] : s_q)
                                    {
                                        sMin = std::min(sMin, s);
                                        sMax = std::max(sMax, s);
                                    }
                                    float L = sMax - sMin;
                                    if(L <= 0.f)
                                        return -1.f;
                                    float Qtotal = 0.f, Qtail = 0.f;
                                    float sCut = sMin + 0.8f * L;
                                    for(auto [s, q] : s_q)
                                    {
                                        Qtotal += q;
                                        if(s >= sCut)
                                            Qtail += q;
                                    }
                                    if(Qtotal <= 0.f)
                                        return -1.f;
                                    return Qtotal / Qtail;
                                },
                                {"MergerData", "TPCData", "MinThetaIdx"});
    ;

    // cout for debuging
    // dfCharge.Foreach([](double th1, double th2)
    //                  { std::cout << "Theta Light 1: " << th1 << " , Theta Light 2: " << th2 << std::endl; },
    //                  {"ThetaLight1", "ThetaLight2"});

    // Debug scaling and definition of distance
    auto dfChargeDiff =
        dfCharge
            .Define("DeltaChargeLight1", [](float ql1, float qb) { return qb - ql1; }, {"ChargeLight1", "ChargeBeam"})
            .Define("DeltaChargeLight2", [](float ql2, float qb) { return qb - ql2; }, {"ChargeLight2", "ChargeBeam"})
            .Define("DeltaChargeLight1aveQfrac", [](float ql1frac, float qb) { return qb - ql1frac; },
                    {"ChargeLight1AveQfrac", "ChargeBeamAveQfrac"})
            .Define("DeltaChargeLight2aveQfrac", [](float ql2frac, float qb) { return qb - ql2frac; },
                    {"ChargeLight2AveQfrac", "ChargeBeamAveQfrac"})
            .Define("DistBeam",
                    [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                    {
                        int beamIdx = m.fBeamIdx;
                        auto clusters = tpc.fClusters;
                        clusters[beamIdx].SortAlongDir(clusters[beamIdx].GetLine().GetDirection());
                        auto voxels = clusters[beamIdx].GetRefToVoxels();
                        ROOT::Math::XYZPointF firstPos = voxels[0].GetPosition();
                        ROOT::Math::XYZPointF lastPos = voxels[voxels.size() - 1].GetPosition();
                        // Scale points to mm
                        ScalePoint(firstPos, 2, 2.84032);
                        ScalePoint(lastPos, 2, 2.84032);
                        auto distance = (firstPos - m.fRP).R();
                        return distance;
                    },
                    {"MergerData", "TPCData"})
            .Define("BeamFirstX",
                    [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                    {
                        int beamIdx = m.fBeamIdx;
                        auto clusters = tpc.fClusters;
                        clusters[beamIdx].SortAlongDir(clusters[beamIdx].GetLine().GetDirection());
                        auto voxels = clusters[beamIdx].GetRefToVoxels();
                        ROOT::Math::XYZPointF firstPos = voxels[0].GetPosition();
                        // Scale points to mm
                        ScalePoint(firstPos, 2, 2.84032);
                        return firstPos.X();
                    },
                    {"MergerData", "TPCData"})
            .Define("BeamLastX",
                    [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                    {
                        int beamIdx = m.fBeamIdx;
                        auto clusters = tpc.fClusters;
                        clusters[beamIdx].SortAlongDir(clusters[beamIdx].GetLine().GetDirection());
                        auto voxels = clusters[beamIdx].GetRefToVoxels();
                        ROOT::Math::XYZPointF lastPos = voxels[voxels.size() - 1].GetPosition();
                        // Scale points to mm
                        ScalePoint(lastPos, 2, 2.84032);
                        return lastPos.X();
                    },
                    {"MergerData", "TPCData"})
            .Define("RPX", [](ActRoot::MergerData& m) { return m.fRP.X(); }, {"MergerData"});

    // Let's try to use the silicon information to try to identify events of the lower Q/L particle
    dfChargeDiff.Foreach(
        [](ActRoot::MergerData& m)
        {
            std::cout << "------------ Event Debug Info -----------" << std::endl;
            std::cout << "Event Silicon Layers hit: " << m.fSilLayers.size() << std::endl;
            for(auto layer : m.fSilLayers)
            {
                std::cout << "Layer: " << layer << std::endl;
            }
            std::cout<<"----------------------------------------" << std::endl;
        },
        {"MergerData"});


    // Check some events if needed
    // std::ofstream out("./Outputs/events_theta_equal.dat");
    // dfChargeDiff.Foreach(
    //     [&](const ActRoot::MergerData& m, double thetaMax, double thetaMin)
    //     {
    //         double tol = 1e-3;
    //         if(std::abs(thetaMax - thetaMin) < tol)
    //         {
    //             m.Stream(out);
    //         }
    //     },
    //     {"MergerData", "ThetaQLengthMax", "ThetaQLengthMin"});
    // out.close();

    // Create histograms to visualize the charge deposition of the three kind of tracks
    // One track on each color
    auto hChargeBeam = dfCharge.Histo1D(
        {"hChargeBeam", "Charge Deposition Beam Particle;Charge/Distance;Counts", 150, 0, 2000}, "ChargeBeam");
    auto hChargeLight1 = dfCharge.Histo1D(
        {"hChargeLight1", "Charge Deposition Light Particle 1;Charge/Distance;Counts", 400, 0, 2000}, "ChargeLight1");
    auto hChargeLight2 = dfCharge.Histo1D(
        {"hChargeLight2", "Charge Deposition Light Particle 2;Charge/Distance;Counts", 400, 0, 2000}, "ChargeLight2");

    auto hDeltaLight1 = dfChargeDiff.Histo1D(
        {"hDeltaLight1", "Charge difference Light1 - Beam;#Delta(Charge/Distance);Counts", 100, -2000, 2000},
        "DeltaChargeLight1");
    auto hDeltaLight2 = dfChargeDiff.Histo1D(
        {"hDeltaLight2", "Charge difference Light2 - Beam;#Delta(Charge/Distance);Counts", 100, -2000, 2000},
        "DeltaChargeLight2");

    auto hDistBeam =
        dfChargeDiff.Histo1D({"hDistBeam", "Distance Beam Particle;Distance [mm];Counts", 400, 0, 200}, "DistBeam");

    auto hBeamFirstX =
        dfChargeDiff.Histo1D({"hBeamFirstX", "Beam First X Position;X [mm];Counts", 400, 0, 200}, "BeamFirstX");
    auto hBeamLastX =
        dfChargeDiff.Histo1D({"hBeamLastX", "Beam Last X Position;X [mm];Counts", 400, 0, 200}, "BeamLastX");
    auto hRPX = dfChargeDiff.Histo1D({"hRPX", "Reaction Point X;X [mm];Counts", 400, 0, 200}, "RPX");

    auto hChargeLight1AveQfrac =
        dfCharge.Histo1D({"hChargeLight1AveQfrac", "Light1 Average Q fraction;Q total / Q tail;Counts", 100, 0, 10},
                         "ChargeLight1AveQfrac");
    auto hChargeLight2AveQfrac =
        dfCharge.Histo1D({"hChargeLight2AveQfrac", "Light2 Average Q fraction;Q total / Q tail;Counts", 100, 0, 10},
                         "ChargeLight2AveQfrac");

    auto hDeltaLight1AveQfrac = dfChargeDiff.Histo1D(
        {"hDeltaLight1AveQfrac", "Average Q fraction difference Light1 - Beam;#Delta(Q total / Q tail);Counts", 200,
         -10, 10},
        "DeltaChargeLight1aveQfrac");
    auto hDeltaLight2AveQfrac = dfChargeDiff.Histo1D(
        {"hDeltaLight2AveQfrac", "Average Q fraction difference Light2 - Beam;#Delta(Q total / Q tail);Counts", 200,
         -10, 10},
        "DeltaChargeLight2aveQfrac");

    auto hThetaLight1 = dfCharge.Histo1D(
        {"hThetaLight1", "Angle Light Particle 1;#Theta Light 1 [deg];Counts", 60, -30, 30}, "ThetaLight1");
    auto hThetaLight2 = dfCharge.Histo1D(
        {"hThetaLight2", "Angle Light Particle 2;#Theta Light 2 [deg];Counts", 60, -30, 30}, "ThetaLight2");
    auto hThetaLight1Light2 = dfCharge.Histo2D(
        {"hThetaLight1Light2", "Angle Light Particles 1 vs 2;#Theta Light 1 [deg];#Theta Light 2 [deg]", 60, -30, 30,
         60, -30, 30},
        "ThetaLight1", "ThetaLight2");

    auto hThetaQLengthMin = dfCharge.Histo1D(
        {"hThetaQLengthMin", "Angle at Min Q/Length;#Theta at Min Q/L [deg];Counts", 60, -30, 30}, "ThetaQLengthMin");
    auto hThetaQLengthMax = dfCharge.Histo1D(
        {"hThetaQLengthMax", "Angle at Max Q/Length;#Theta at Max Q/L [deg];Counts", 60, -30, 30}, "ThetaQLengthMax");
    auto hThetaQLengthMinVsMax = dfCharge.Histo2D(
        {"hThetaQLengthMinVsMax", "Angle at Min vs Max Q/Length;#Theta at Min Q/L [deg];#Theta at Max Q/L [deg]", 60,
         -30, 30, 60, -30, 30},
        "ThetaQLengthMin", "ThetaQLengthMax");

    auto hQLengthMax = dfCharge.Histo1D(
        {"hQLengthMax", "Max Q/Length;Max Q/Length [Charge/Distance];Counts", 150, 0, 2000}, "QLengthMax");
    auto hQLengthMin = dfCharge.Histo1D(
        {"hQLengthMin", "Min Q/Length;Min Q/Length [Charge/Distance];Counts", 150, 0, 2000}, "QLengthMin");
    auto hQLengthMinVsMax = dfCharge.Histo2D(
        {"hQLengthMinVsMax", "Min vs Max Q/Length;Min Q/Length [Charge/Distance];Max Q/Length [Charge/Distance]", 150,
         0, 2000, 150, 0, 2000},
        "QLengthMin", "QLengthMax");

    // Plot
    auto c1 {new TCanvas("c1", "Charge Deposition", 800, 600)};
    c1->cd();
    hChargeBeam->SetLineColor(kRed);
    hChargeLight1->SetLineColor(kBlue);
    hChargeLight2->SetLineColor(kGreen);
    hChargeBeam->DrawClone();
    hChargeLight1->DrawClone("SAME");
    hChargeLight2->DrawClone("SAME");
    // c1.BuildLegend();
    auto c2 {new TCanvas("c2", "Delta Charge Deposition", 800, 600)};
    c2->cd();
    hDeltaLight1->SetLineColor(kBlue);
    hDeltaLight2->SetLineColor(kGreen);
    hDeltaLight1->DrawClone();
    hDeltaLight2->DrawClone("SAME");
    //
    auto c3 {new TCanvas("c3", "Distance Beam Particle", 800, 600)};
    c3->cd();
    hDistBeam->DrawClone();
    hRPX->SetLineColor(kRed);
    hRPX->DrawClone("same");
    //
    auto c4 {new TCanvas("c4", "Beam X Positions", 800, 600)};
    c4->cd();
    hBeamFirstX->SetLineColor(kRed);
    hBeamLastX->SetLineColor(kBlue);
    hBeamFirstX->DrawClone();
    hBeamLastX->DrawClone("SAME");
    //
    auto c5 {new TCanvas("c5", "Light1 Average Q fraction", 800, 600)};
    c5->cd();
    hChargeLight1AveQfrac->DrawClone();
    hChargeLight2AveQfrac->SetLineColor(kRed);
    hChargeLight2AveQfrac->DrawClone("SAME");
    //
    auto c6 {new TCanvas("c6", "Delta Light Average Q fraction", 800, 600)};
    c6->cd();
    hDeltaLight1AveQfrac->DrawClone();
    hDeltaLight2AveQfrac->SetLineColor(kRed);
    hDeltaLight2AveQfrac->DrawClone("SAME");
    //
    auto c7 {new TCanvas("c7", "Light Particles Angles", 800, 600)};
    c7->Divide(3, 1);
    c7->cd(1);
    hThetaLight1->DrawClone();
    c7->cd(2);
    hThetaLight2->SetLineColor(kRed);
    hThetaLight2->DrawClone("SAME");
    c7->cd(3);
    hThetaLight1Light2->DrawClone();
    //
    auto c8 {new TCanvas("c8", "Theta at Min and Max Q/Length", 800, 600)};
    c8->Divide(3, 1);
    c8->cd(1);
    hThetaQLengthMin->DrawClone();
    c8->cd(2);
    hThetaQLengthMax->SetLineColor(kRed);
    hThetaQLengthMax->DrawClone("SAME");
    c8->cd(3);
    hThetaQLengthMinVsMax->DrawClone();
    //
    auto c9 {new TCanvas("c9", "Min and Max Q/Length", 800, 600)};
    c9->Divide(2, 1);
    c9->cd(1);
    hQLengthMin->DrawClone();
    hQLengthMax->SetLineColor(kRed);
    hQLengthMax->DrawClone("SAME");
    hChargeBeam->SetLineColor(kGreen + 1);
    hChargeBeam->DrawClone("SAME");
    c9->cd(2);
    hQLengthMinVsMax->DrawClone();
}