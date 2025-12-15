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

// PASA ALGO COA DISTANCIA, ANTES SUMANDO A DISTANCIA ENTRE CADA VOXEL DABA UNHA ACUMULACION DE CARGA PARA O BEAM , EPRO
// AGORA NO PENSO QUE TAMPOUCO ESTÁ ESCALADA A MAGNITUDES FISICAS

// float distance = 0.0f;
// ROOT::Math::XYZPointF prevPos = voxels[0].GetPosition();
// // Loop over voxels to calculate total charge and distance
// for(const auto& voxel : voxels)
// {
//     totalCharge += voxel.GetCharge();
//     auto currentPos = voxel.GetPosition();
//     distance += std::sqrt(std::pow(currentPos.X() - prevPos.X(), 2) + std::pow(currentPos.Y() - prevPos.Y(), 2) +
//                           std::pow(currentPos.Z() - prevPos.Z(), 2));
//     prevPos = currentPos;
// }

void ScalePoint(ROOT::Math::XYZPointF& point, float xy, float z, bool addOffset = false)
{
    if(addOffset) // when converting a bin point to physical units which wasnt already corrected
        point += ROOT::Math::XYZVector {0.5, 0.5, 0.5};
    point.SetX(point.X() * xy);
    point.SetY(point.Y() * xy);
    point.SetZ(point.Z() * z);
}

void CheckChargeDeposition()
{
    // ROOT::DisableImplicitMT();
    ROOT::RDataFrame df {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_7Li_d_p_filtered.root"};

    // Start with only silicon events
    auto dfFilter = df.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == true; }, {"MergerData"});

    //  Define the charge / dist for each particle except light particle detected in silicon
    auto dfCharge = dfFilter
                        .Define("ChargeBeam",
                                [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                                {
                                    int beamIdx = m.fBeamIdx;
                                    auto clusters = tpc.fClusters;
                                    clusters[beamIdx].SortAlongDir(clusters[beamIdx].GetLine().GetDirection());
                                    auto voxels = clusters[beamIdx].GetRefToVoxels();
                                    // Define variables for total charge and distance
                                    float totalCharge = 0.0f;
                                    ROOT::Math::XYZPointF firstPos = voxels[0].GetPosition();
                                    ROOT::Math::XYZPointF lastPos = voxels[voxels.size() - 1].GetPosition();
                                    // Scale points to mm
                                    ScalePoint(firstPos, 2, 2.84032);
                                    ScalePoint(lastPos, 2, 2.84032);
                                    auto distance = (lastPos - firstPos).R();
                                    // Loop over voxels to calculate total charge and distance
                                    for(const auto& voxel : voxels)
                                    {
                                        totalCharge += voxel.GetCharge();
                                    }
                                    // Avoid division by zero
                                    if(distance == 0.0f)
                                        return 0.0f;
                                    return totalCharge / distance; // Charge per unit distance
                                },
                                {"MergerData", "TPCData"})
                        .Define("ChargeLight1",
                                [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                                {
                                    int lightIdx = m.fLightIdx;
                                    int beamIdx = m.fBeamIdx;
                                    auto clusters = tpc.fClusters;
                                    // Get index diferent than lightIdx and beamIdx to analyse first diferent that has
                                    // higher angle respect to beam
                                    int idx = -1;
                                    double maxAngle = -1.0;
                                    for(int i = 0; i < clusters.size(); i++)
                                    {
                                        if(i != lightIdx && i != beamIdx)
                                        {
                                            // Calculate the angle between the beam and the current cluster
                                            auto beamDirection = clusters[beamIdx].GetLine().GetDirection();
                                            auto currentDirection = clusters[i].GetLine().GetDirection();
                                            double angle = beamDirection.Dot(currentDirection);
                                            if(angle > maxAngle)
                                            {
                                                maxAngle = angle;
                                                idx = i;
                                            }
                                        }
                                    }
                                    clusters[idx].SortAlongDir(clusters[idx].GetLine().GetDirection());
                                    auto voxels = clusters[idx].GetRefToVoxels();
                                    // Define variables for total charge and distance
                                    float totalCharge = 0.0f;
                                    ROOT::Math::XYZPointF firstPos = voxels[0].GetPosition();
                                    ROOT::Math::XYZPointF lastPos = voxels[voxels.size() - 1].GetPosition();
                                    // Scale points to mm
                                    ScalePoint(firstPos, 2, 2.84032);
                                    ScalePoint(lastPos, 2, 2.84032);
                                    auto distance = (lastPos - firstPos).R();
                                    // Loop over voxels to calculate total charge
                                    for(const auto& voxel : voxels)
                                    {
                                        totalCharge += voxel.GetCharge();
                                    }
                                    // Avoid division by zero
                                    if(distance == 0.0f)
                                        return 50.0f;
                                    return totalCharge / distance; // Charge per unit distance
                                },
                                {"MergerData", "TPCData"})
                        .Define("ChargeLight2",
                                [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                                {
                                    int lightIdx = m.fLightIdx;
                                    int beamIdx = m.fBeamIdx;
                                    auto clusters = tpc.fClusters;
                                    // Get index diferent than lightIdx and beamIdx to analyse first diferent that has
                                    // lower angle respect to beam
                                    int idx = -1;
                                    double minAngle = 100000;
                                    for(int i = 0; i < clusters.size(); i++)
                                    {
                                        if(i != lightIdx && i != beamIdx)
                                        {
                                            // Calculate the angle between the beam and the current cluster
                                            auto beamDirection = clusters[beamIdx].GetLine().GetDirection();
                                            auto currentDirection = clusters[i].GetLine().GetDirection();
                                            double angle = beamDirection.Dot(currentDirection);
                                            if(angle < minAngle)
                                            {
                                                minAngle = angle;
                                                idx = i;
                                            }
                                        }
                                    }
                                    clusters[idx].SortAlongDir(clusters[idx].GetLine().GetDirection());
                                    auto voxels = clusters[idx].GetRefToVoxels();
                                    // Define variables for total charge and distance
                                    float totalCharge = 0.0f;
                                    ROOT::Math::XYZPointF firstPos = voxels[0].GetPosition();
                                    ROOT::Math::XYZPointF lastPos = voxels[voxels.size() - 1].GetPosition();
                                    // Scale points to mm
                                    ScalePoint(firstPos, 2, 2.84032);
                                    ScalePoint(lastPos, 2, 2.84032);
                                    auto distance = (lastPos - firstPos).R();
                                    // Loop over voxels to calculate total charge and distance
                                    for(const auto& voxel : voxels)
                                    {
                                        totalCharge += voxel.GetCharge();
                                    }
                                    // Avoid division by zero
                                    if(distance == 0.0f)
                                        return 50.0f;
                                    if(totalCharge == 0.0f)
                                        return 50.0f;
                                    return totalCharge / distance; // Charge per unit distance
                                },
                                {"MergerData", "TPCData"})
                        .Define("ChargeBeamAveQfrac",
                                [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                                {
                                    int beamIdx = m.fBeamIdx;
                                    auto clusters = tpc.fClusters;

                                    auto& cluster = clusters[beamIdx];
                                    auto voxels = cluster.GetRefToVoxels();

                                    if(voxels.size() < 5)
                                        return -100.f;

                                    // --- dirección longitudinal del beam
                                    auto dir = cluster.GetLine().GetDirection().Unit();

                                    // --- origen longitudinal
                                    ROOT::Math::XYZPointF r0 = voxels.front().GetPosition();

                                    // --- coordenada longitudinal proyectada de cada voxel
                                    std::vector<std::pair<float, float>> s_q;
                                    s_q.reserve(voxels.size());

                                    for(const auto& v : voxels)
                                    {
                                        auto r = v.GetPosition();
                                        float s = (r - r0).Dot(dir);
                                        s_q.emplace_back(s, v.GetCharge());
                                    }

                                    // --- rango longitudinal
                                    float sMin = 1e9f;
                                    float sMax = -1e9f;

                                    for(const auto& [s, q] : s_q)
                                    {
                                        sMin = std::min(sMin, s);
                                        sMax = std::max(sMax, s);
                                    }

                                    float L = sMax - sMin;
                                    if(L <= 0.f)
                                        return -100.f;

                                    // --- integrar carga total y carga en el último 30 %
                                    float Qtotal = 0.f;
                                    float Qtail = 0.f;

                                    float sCut = sMin + 0.7f * L;

                                    for(const auto& [s, q] : s_q)
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
                                [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                                {
                                    int lightIdx = m.fLightIdx;
                                    int beamIdx = m.fBeamIdx;
                                    auto clusters = tpc.fClusters;

                                    // --- seleccionar el cluster light menos colineal con el beam
                                    int idx = -1;
                                    double maxAngle = -1.0;
                                    for(int i = 0; i < clusters.size(); i++)
                                    {
                                        if(i != lightIdx && i != beamIdx)
                                        {
                                            // Calculate the angle between the beam and the current cluster
                                            auto beamDirection = clusters[beamIdx].GetLine().GetDirection();
                                            auto currentDirection = clusters[i].GetLine().GetDirection();
                                            double angle = beamDirection.Dot(currentDirection);
                                            if(angle > maxAngle)
                                            {
                                                maxAngle = angle;
                                                idx = i;
                                            }
                                        }
                                    }

                                    if(idx < 0)
                                        return -1.f;

                                    auto& cluster = clusters[idx];
                                    cluster.SortAlongDir(cluster.GetLine().GetDirection());
                                    auto voxels = cluster.GetRefToVoxels();

                                    if(voxels.size() < 5)
                                        return -1.f;

                                    // --- dirección longitudinal del track
                                    auto dir = cluster.GetLine().GetDirection().Unit();

                                    // --- origen longitudinal (primer voxel)
                                    ROOT::Math::XYZPointF r0 = voxels.front().GetPosition();

                                    // --- calcular coordenada longitudinal proyectada s
                                    std::vector<std::pair<float, float>> s_q;
                                    s_q.reserve(voxels.size());

                                    for(const auto& v : voxels)
                                    {
                                        auto r = v.GetPosition();
                                        float s = (r - r0).Dot(dir);
                                        s_q.emplace_back(s, v.GetCharge());
                                    }

                                    // --- rango longitudinal
                                    float sMin = 1e9f;
                                    float sMax = -1e9f;

                                    for(const auto& [s, q] : s_q)
                                    {
                                        sMin = std::min(sMin, s);
                                        sMax = std::max(sMax, s);
                                    }

                                    float L = sMax - sMin;
                                    if(L <= 0.f)
                                        return -1.f;

                                    // --- integrar carga total y carga en el último 30 %
                                    float Qtotal = 0.f;
                                    float Qtail = 0.f;

                                    float sCut = sMin + 0.8f * L;

                                    for(const auto& [s, q] : s_q)
                                    {
                                        Qtotal += q;
                                        if(s >= sCut)
                                            Qtail += q;
                                    }

                                    if(Qtotal <= 0.f)
                                        return -1.f;

                                    return Qtotal / Qtail;
                                },
                                {"MergerData", "TPCData"})
                        .Define("ChargeLight2AveQfrac",
                                [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                                {
                                    int lightIdx = m.fLightIdx;
                                    int beamIdx = m.fBeamIdx;
                                    auto clusters = tpc.fClusters;

                                    // --- seleccionar el cluster light más colineal con el beam
                                    int idx = -1;
                                    double minAngle = 100000;
                                    for(int i = 0; i < clusters.size(); i++)
                                    {
                                        if(i != lightIdx && i != beamIdx)
                                        {
                                            // Calculate the angle between the beam and the current cluster
                                            auto beamDirection = clusters[beamIdx].GetLine().GetDirection();
                                            auto currentDirection = clusters[i].GetLine().GetDirection();
                                            double angle = beamDirection.Dot(currentDirection);
                                            if(angle < minAngle)
                                            {
                                                minAngle = angle;
                                                idx = i;
                                            }
                                        }
                                    }

                                    if(idx < 0)
                                        return -1.f;

                                    auto& cluster = clusters[idx];
                                    cluster.SortAlongDir(cluster.GetLine().GetDirection());
                                    auto voxels = cluster.GetRefToVoxels();

                                    if(voxels.size() < 5)
                                        return -1.f;

                                    // --- dirección longitudinal del track
                                    auto dir = cluster.GetLine().GetDirection().Unit();

                                    // --- origen longitudinal (primer voxel)
                                    ROOT::Math::XYZPointF r0 = voxels.front().GetPosition();

                                    // --- calcular coordenada longitudinal proyectada s
                                    std::vector<std::pair<float, float>> s_q;
                                    s_q.reserve(voxels.size());

                                    for(const auto& v : voxels)
                                    {
                                        auto r = v.GetPosition();
                                        float s = (r - r0).Dot(dir);
                                        s_q.emplace_back(s, v.GetCharge());
                                    }

                                    // --- rango longitudinal
                                    float sMin = 1e9f;
                                    float sMax = -1e9f;

                                    for(const auto& [s, q] : s_q)
                                    {
                                        sMin = std::min(sMin, s);
                                        sMax = std::max(sMax, s);
                                    }

                                    float L = sMax - sMin;
                                    if(L <= 0.f)
                                        return -1.f;

                                    // --- integrar carga total y carga en el último 30 %
                                    float Qtotal = 0.f;
                                    float Qtail = 0.f;

                                    float sCut = sMin + 0.8f * L;

                                    for(const auto& [s, q] : s_q)
                                    {
                                        Qtotal += q;
                                        if(s >= sCut)
                                            Qtail += q;
                                    }

                                    if(Qtotal <= 0.f)
                                        return -1.f;

                                    return Qtotal / Qtail;
                                },
                                {"MergerData", "TPCData"});

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

    // Let's try to use the silicon information to try to identify some events
    

    // Check some events if needed
    // std::ofstream out("./Outputs/rpx_greater_120.dat");
    // dfChargeDiff.Foreach(
    //     [&](ActRoot::MergerData& m, float rpx)
    //     {
    //         if(rpx > 120)
    //             m.Stream(out);
    //     },
    //     {"MergerData", "RPX"});
    // out.close();

    // Create histograms to visualize the charge deposition of the three kind of tracks
    // One track on each color
    auto hChargeBeam = dfCharge.Histo1D(
        {"hChargeBeam", "Charge Deposition Beam Particle;Charge/Distance;Counts", 400, 0, 2000}, "ChargeBeam");
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
        dfCharge.Histo1D({"hChargeLight1AveQfrac", "Light1 Average Q fraction;Q total / Q tail;Counts", 100, 0, 1},
                         "ChargeLight1AveQfrac");
    auto hChargeLight2AveQfrac =
        dfCharge.Histo1D({"hChargeLight2AveQfrac", "Light2 Average Q fraction;Q total / Q tail;Counts", 100, 0, 1},
                         "ChargeLight2AveQfrac");

    auto hDeltaLight1AveQfrac = dfChargeDiff.Histo1D(
        {"hDeltaLight1AveQfrac", "Average Q fraction difference Light1 - Beam;#Delta(Q total / Q tail);Counts", 200, -10,
         10},
        "DeltaChargeLight1aveQfrac");
    auto hDeltaLight2AveQfrac = dfChargeDiff.Histo1D(
        {"hDeltaLight2AveQfrac", "Average Q fraction difference Light2 - Beam;#Delta(Q total / Q tail);Counts", 200, -10,
         10},
        "DeltaChargeLight2aveQfrac");

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
}