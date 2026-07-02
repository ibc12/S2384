#include "ActContinuity.h"
#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActTPCData.h"
#include "ActTPCParameters.h"
#include "ActVoxel.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TH2.h"
#include "TMath.h"

#include "Math/Point3Dfwd.h"

#include <cstdint>
#include <fstream>
#include <iostream>
#include <unordered_set>

// COnvert position into a int for a fast lookup
static inline int64_t HashCell(int ix, int iy, int iz)
{
    return ((int64_t)ix << 22) | ((int64_t)iy << 11) | (int64_t)iz;
}
std::vector<ROOT::Math::XYZPoint>
FindSandwichedHoles(const std::vector<ActRoot::Voxel>& voxels, double dx, double dy, double dz, int minNeighbors = 26)
{
    std::vector<ROOT::Math::XYZPoint> holes;

    std::unordered_set<int64_t> occupied;
    occupied.reserve(voxels.size() * 2);
    for(auto& voxel : voxels)
    {
        auto pos = voxel.GetPosition();
        occupied.insert(HashCell(int(pos.X()), int(pos.Y()), int(pos.Z())));
    }

    // los 26 vecinos de una celda (todas las combinaciones -1,0,1 excepto (0,0,0))
    std::vector<std::array<int, 3>> neighbors26;
    for(int dx_ = -1; dx_ <= 1; dx_++)
        for(int dy_ = -1; dy_ <= 1; dy_++)
            for(int dz_ = -1; dz_ <= 1; dz_++)
            {
                if(dx_ == 0 && dy_ == 0 && dz_ == 0)
                    continue;
                neighbors26.push_back({dx_, dy_, dz_});
            }

    std::unordered_set<int64_t> candidates; // celdas vacías a examinar
    std::unordered_set<int64_t> alreadyFlagged;

    // Paso 1: reunir todas las celdas vacías que son vecinas de al menos un voxel ocupado
    for(auto& voxel : voxels)
    {
        auto pos = voxel.GetPosition();
        auto [ix, iy, iz] = std::make_tuple(int(pos.X()), int(pos.Y()), int(pos.Z()));

        for(auto& d : neighbors26)
        {
            int cx = ix + d[0], cy = iy + d[1], cz = iz + d[2];
            int64_t candKey = HashCell(cx, cy, cz);
            if(occupied.count(candKey))
                continue; // no es un hueco, ya está ocupada
            candidates.insert(candKey);
        }
    }

    // Paso 2: para cada celda vacía candidata, contar cuántos de sus 26 vecinos están ocupados
    for(auto& candKey : candidates)
    {
        if(alreadyFlagged.count(candKey))
            continue;

        // deshacer el hash para recuperar ix, iy, iz
        // (ajusta esto a los shifts reales que estés usando en HashCell)
        int cz = candKey & 0x7FF;         // últimos 11 bits
        int cy = (candKey >> 11) & 0x7FF; // siguientes 11 bits
        int cx = (candKey >> 22) & 0x7FF; // últimos 11 bits restantes

        int occupiedCount = 0;
        for(auto& d : neighbors26)
        {
            int64_t nKey = HashCell(cx + d[0], cy + d[1], cz + d[2]);
            if(occupied.count(nKey))
                occupiedCount++;
        }

        if(occupiedCount < minNeighbors)
            continue;

        alreadyFlagged.insert(candKey);
        ROOT::Math::XYZPoint h(cx * dx, cy * dy, cz * dz);
        holes.push_back(h);
    }

    return holes;
}

void TracksMissingPadsZNeighbours()
{
    // ROOT::EnableImplicitMT();

    // ROOT::RDataFrame df {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_7Li_d_d_filtered.root"};
    // ROOT::RDataFrame df {"PreProcessed_Tree", "../../PostAnalysis/Outputs/tree_preprocess_F_11Li.root"};

    ActRoot::DataManager dataman {"../../configs/data_7Li.conf", ActRoot::ModeType::EMerge};
    dataman.SetRuns(69, 69);
    auto chain {dataman.GetChain()};
    // Add friends if necessary
    auto friend1 {dataman.GetChain(ActRoot::ModeType::EFilter)};
    chain->AddFriend(friend1.get());
    auto friend2 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(friend2.get());

    ROOT::RDataFrame df {*chain};
    // auto def {
    //     df.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"})}; // only L1
    // auto def {df.Filter([](ActRoot::ModularData& m) { return m.Get("GATCONF") == 8; }, {"ModularData"})}; // only L1
    // Filter to get just one event
    auto def {df.Filter([](ActRoot::MergerData& m) { return m.fEntry == 5876; }, {"MergerData"})}; // only L1 (e5889)
    // Let's get a goos L1 event with high phi, and compare it with a missing voxels event
    // We will focus on one event that has a hole, to see if we identify it.
    auto defZ = def.Define("Z",
                           [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                           {
                               auto lightIdx = m.fLightIdx;
                               auto lightCl = tpc.fClusters[lightIdx];
                               auto voxels = lightCl.GetVoxels();
                               std::vector<float> zPositions;
                               for(auto& v : voxels)
                               {
                                   zPositions.push_back(v.GetPosition().Z());
                               }
                               return zPositions;
                           },
                           {"MergerData", "TPCData"})
                    .Define("nVoxels",
                            [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                            {
                                auto lightIdx = m.fLightIdx;
                                auto lightCl = tpc.fClusters[lightIdx];
                                return lightCl.GetVoxels().size();
                            },
                            {"MergerData", "TPCData"});
    auto hPosZ = defZ.Histo1D({"hPosZ", "Z positions of voxels;Z [tb];Counts", 512, 0, 512}, "Z");
    auto hNVoxels =
        defZ.Histo1D({"hNVoxels", "Number of voxels in light cluster;# Voxels;Counts", 100, 0, 1000}, "nVoxels");
    auto* c = new TCanvas("c", "Z positions of voxels", 800, 600);
    c->Divide(2, 1);
    c->cd(1);
    hNVoxels->DrawClone();
    c->cd(2);
    hPosZ->DrawClone();

    auto defHoles = defZ.Define("holes",
                                [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                                {
                                    auto lightIdx = m.fLightIdx;
                                    auto lightCl = tpc.fClusters[lightIdx];
                                    auto voxels = lightCl.GetVoxels();
                                    // dx, dy, dz: your voxel pitch (pad pitch in x/y, time-bucket size in z)
                                    double dx = 1.0; // mm
                                    double dy = 1.0; // mm
                                    double dz = 1.0; // 4tb = 1btb
                                    return FindSandwichedHoles(voxels, dx, dy, dz, 8);
                                },
                                {"MergerData", "TPCData"})
                        .Define("holesZ",
                                [](const std::vector<ROOT::Math::XYZPoint>& holes)
                                {
                                    std::vector<float> zPositions;
                                    for(auto& h : holes)
                                    {
                                        zPositions.push_back(h.Z());
                                    }
                                    return zPositions;
                                },
                                {"holes"});
    auto hHolesZ = defHoles.Histo1D({"hHolesZ", "Z positions of holes;Z [tb];Counts", 128, 0, 128}, "holesZ");
    auto* c2 = new TCanvas("c2", "Z positions of holes", 800, 600);
    hHolesZ->DrawClone();
    // dx, dy, dz: your voxel pitch (pad pitch in x/y, time-bucket size in z)
    // Used to quantize positions onto an integer grid so we can do exact
    // neighbor lookups instead of distance comparisons.

    // Plot the xy xz progjections of the voxels of the light cluster
    auto hXY = new TH2F("hXY", "XY projection of voxels;X [mm];Y [mm]", 128, 0, 128, 128, 0, 128);
    auto hXZ = new TH2F("hXZ", "XZ projection of voxels;X [mm];Z [tb]", 128, 0, 128, 128, 0, 128);
    auto hYZ = new TH2F("hYZ", "YZ projection of voxels;Y [mm];Z [tb]", 128, 0, 128, 128, 0, 128);

    std::vector<double> holesX, holesY, holesZv; // holesZv para no chocar con la var "holesZ" de arriba

    defHoles.Foreach(
        [&hXY, &hXZ, &hYZ, &holesX, &holesY, &holesZv](ActRoot::MergerData& m, ActRoot::TPCData& tpc,
                                                       const std::vector<ROOT::Math::XYZPoint>& holes)
        {
            auto lightIdx = m.fLightIdx;
            std::cout << "Light index: " << lightIdx << std::endl;
            auto lightCl = tpc.fClusters[lightIdx];
            auto voxels = lightCl.GetVoxels();
            for(auto& v : voxels)
            {
                auto pos = v.GetPosition();
                std::cout << "Voxel position: (" << pos.X() << ", " << pos.Y() << ", " << pos.Z() << ")" << std::endl;
                hXY->Fill(pos.X(), pos.Y());
                hXZ->Fill(pos.X(), pos.Z());
                hYZ->Fill(pos.Y(), pos.Z());
            }
            // Acumular los holes de este evento
            for(auto& h : holes)
            {
                std::cout << "Hole position: (" << h.X() << ", " << h.Y() << ", " << h.Z() << ")" << std::endl;
                holesX.push_back(h.X() + 0.5); // +0.5 para centrar en el bin
                holesY.push_back(h.Y() + 0.5);
                holesZv.push_back(h.Z() + 0.5);
            }
        },
        {"MergerData", "TPCData", "holes"});
    // Construir los TGraph de holes para cada proyección
    auto* gXY = new TGraph(holesX.size(), holesX.data(), holesY.data());
    auto* gXZ = new TGraph(holesX.size(), holesX.data(), holesZv.data());
    auto* gYZ = new TGraph(holesY.size(), holesY.data(), holesZv.data());

    for(auto* g : {gXY, gXZ, gYZ})
    {
        g->SetMarkerStyle(29); // estrella
        g->SetMarkerSize(1.5);
        g->SetMarkerColor(kRed);
    }

    auto* c3 = new TCanvas("c3", "Projections of voxels + holes", 1200, 400);
    c3->Divide(3, 1);
    c3->cd(1);
    hXY->DrawClone("colz");
    gXY->Draw("P same");
    c3->cd(2);
    hXZ->DrawClone("colz");
    gXZ->Draw("P same");
    c3->cd(3);
    hYZ->DrawClone("colz");
    gYZ->Draw("P same");
}