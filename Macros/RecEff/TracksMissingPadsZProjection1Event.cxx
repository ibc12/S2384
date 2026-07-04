#include "ActContinuity.h"
#include "ActDataManager.h"
#include "ActLine.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActTPCData.h"
#include "ActTPCParameters.h"
#include "ActVoxel.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"

#include "Math/Point3Dfwd.h"
#include "Math/Vector3D.h"

#include <algorithm>
#include <iostream>
#include <vector>

// Ahora devolvemos, para cada voxel, tanto su coordenada proyectada s
// como su posición 3D real, ordenados por s. Así podemos luego recuperar
// las posiciones reales de los voxeles que rodean el gap más grande.
struct ProjectedPoint
{
    double s;
    ROOT::Math::XYZPointF pos3D;
};

std::vector<ProjectedPoint>
GetProjectedPositions(ActRoot::Cluster lightCl, const ROOT::Math::XYZPointF& ref, double driftFactor = 1.)
{
    // auto voxels = lightCl.GetRefToVoxels();
    lightCl.ScaleVoxels(2, driftFactor); // This function already puts 0.5 offset
    lightCl.ReFit();
    auto voxels = lightCl.GetRefToVoxels();
    auto line = lightCl.GetRefToLine();
    auto dir = line.GetDirection().Unit();
    lightCl.SortAlongDir(dir); // ordenamos los voxeles a lo largo de la dirección del cluster
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

void TracksMissingPadsZProjection1Event()
{
    ActRoot::DataManager dataman {"../../configs/data_7Li.conf", ActRoot::ModeType::EMerge};
    dataman.SetRuns(69, 69);
    auto chain {dataman.GetChain()};
    auto friend1 {dataman.GetChain(ActRoot::ModeType::EFilter)};
    chain->AddFriend(friend1.get());
    auto friend2 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(friend2.get());

    ActRoot::InputParser parser {};
    parser.ReadFile("../../configs/detector.conf");
    auto driftBlock = parser.GetBlock("Merger");
    auto driftFactor = driftBlock->GetDouble("DriftFactor");

    ROOT::RDataFrame df {*chain};
    auto def {df.Filter([](ActRoot::MergerData& m) { return m.fEntry == 5889; }, {"MergerData"})};

    auto defProj = def.Define("projPoints",
                              [&driftFactor](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                              {
                                  auto lightIdx = m.fLightIdx;
                                  auto lightCl = tpc.fClusters[lightIdx];
                                  auto rp = m.fRP;
                                  return GetProjectedPositions(lightCl, rp, driftFactor);
                              },
                              {"MergerData", "TPCData"});

    auto projVecs = defProj.Take<std::vector<ProjectedPoint>>("projPoints");
    if(projVecs->empty())
    {
        std::cerr << "No events passed the filter!" << std::endl;
        return;
    }
    auto& points = projVecs->at(0); // solo un evento

    std::cout << "Number of projected points: " << points.size() << std::endl;
    std::cout << "s range: [" << points.front().s << ", " << points.back().s << "]" << std::endl;

    // Buscar el gap más grande y guardarnos el índice donde ocurre
    double maxGap = -1;
    size_t maxGapIdx = 0; // el gap está entre points[maxGapIdx-1] y points[maxGapIdx]
    std::cout << "\nConsecutive gaps (s[i] - s[i-1]):" << std::endl;
    for(size_t i = 1; i < points.size(); ++i)
    {
        double gap = points[i].s - points[i - 1].s;
        std::cout << "  i=" << i << "  gap=" << gap << std::endl;
        if(gap > maxGap)
        {
            maxGap = gap;
            maxGapIdx = i;
        }
    }

    ROOT::Math::XYZPointF pBefore, pAfter;
    bool haveMaxGap = (points.size() >= 2);
    if(haveMaxGap)
    {
        pBefore = points[maxGapIdx - 1].pos3D;
        pAfter = points[maxGapIdx].pos3D;
        std::cout << "\nBiggest gap = " << maxGap << " between index " << (maxGapIdx - 1) << " and " << maxGapIdx
                  << std::endl;
        std::cout << "  Before: (" << pBefore.X() << ", " << pBefore.Y() << ", " << pBefore.Z() << ")" << std::endl;
        std::cout << "  After:  (" << pAfter.X() << ", " << pAfter.Y() << ", " << pAfter.Z() << ")" << std::endl;
    }

    // --- Histogramas 1D de s con distintos binnings (igual que antes) ---
    std::vector<double> sPos;
    sPos.reserve(points.size());
    for(auto& p : points)
        sPos.push_back(p.s);

    double sMin = sPos.front();
    double sMax = sPos.back();
    double range = sMax - sMin;
    double pad = 0.05 * range;
    sMin -= pad;
    sMax += pad;

    std::vector<int> nBinsOptions = {50, 100, 200, 400};
    auto* c1 = new TCanvas("c1", "Projected positions - different binnings", 1200, 800);
    c1->Divide(2, 2);
    std::vector<TH1F*> hists;
    for(size_t i = 0; i < nBinsOptions.size(); ++i)
    {
        int nBins = nBinsOptions[i];
        auto* h =
            new TH1F(Form("hProj_%d", nBins), Form("Projected s;s [mm];Counts (bins=%d)", nBins), nBins, sMin, sMax);
        for(auto& s : sPos)
            h->Fill(s);
        hists.push_back(h);
        c1->cd(i + 1);
        h->SetLineWidth(2);
        h->DrawClone("hist");
    }

    std::vector<double> idx(sPos.size());
    for(size_t i = 0; i < sPos.size(); ++i)
        idx[i] = (double)i;

    auto* gProj = new TGraph(sPos.size(), idx.data(), sPos.data());
    gProj->SetTitle("Projected s vs sorted index;index;s [mm]");
    gProj->SetMarkerStyle(20);
    gProj->SetMarkerSize(0.8);
    auto* c2 = new TCanvas("c2", "Projected s vs index (gap visualization)", 800, 600);
    gProj->Draw("AP");

    std::vector<double> gapIdx, gapVal;
    for(size_t i = 1; i < sPos.size(); ++i)
    {
        gapIdx.push_back((double)i);
        gapVal.push_back(sPos[i] - sPos[i - 1]);
    }
    auto* gGaps = new TGraph(gapIdx.size(), gapIdx.data(), gapVal.data());
    gGaps->SetTitle("Consecutive gap size vs index;index;gap = s[i]-s[i-1]  [mm]");
    gGaps->SetMarkerStyle(20);
    gGaps->SetMarkerColor(kRed);
    gGaps->SetMarkerSize(0.8);
    auto* c3 = new TCanvas("c3", "Gap sizes", 800, 600);
    gGaps->Draw("AP");

    // --- Proyecciones XY, XZ, YZ de todos los voxeles, con el gap más grande marcado ---
    auto hXY = new TH2F("hXY", "XY projection;X [mm];Y [mm]", 128, 0, 256, 128, 0, 256);
    auto hXZ = new TH2F("hXZ", "XZ projection;X [mm];Z [tb]", 128, 0, 256, 256, 0, 512);
    auto hYZ = new TH2F("hYZ", "YZ projection;Y [mm];Z [tb]", 128, 0, 256, 128, 0, 512);

    for(auto& p : points)
    {
        hXY->Fill(p.pos3D.X(), p.pos3D.Y());
        hXZ->Fill(p.pos3D.X(), p.pos3D.Z());
        hYZ->Fill(p.pos3D.Y(), p.pos3D.Z());
    }

    // Marcadores para los dos puntos que rodean el gap más grande
    TGraph *gGapXY = nullptr, *gGapXZ = nullptr, *gGapYZ = nullptr;
    if(haveMaxGap)
    {
        double gapXY_x[2] = {pBefore.X(), pAfter.X()};
        double gapXY_y[2] = {pBefore.Y(), pAfter.Y()};
        gGapXY = new TGraph(2, gapXY_x, gapXY_y);

        double gapXZ_x[2] = {pBefore.X(), pAfter.X()};
        double gapXZ_z[2] = {pBefore.Z(), pAfter.Z()};
        gGapXZ = new TGraph(2, gapXZ_x, gapXZ_z);

        double gapYZ_y[2] = {pBefore.Y(), pAfter.Y()};
        double gapYZ_z[2] = {pBefore.Z(), pAfter.Z()};
        gGapYZ = new TGraph(2, gapYZ_y, gapYZ_z);

        for(auto* g : {gGapXY, gGapXZ, gGapYZ})
        {
            g->SetMarkerStyle(29); // estrella
            g->SetMarkerSize(2.0);
            g->SetMarkerColor(kRed);
        }
    }

    auto* c4 = new TCanvas("c4", "Projections + biggest gap", 1200, 400);
    c4->Divide(3, 1);
    c4->cd(1);
    hXY->DrawClone("colz");
    if(gGapXY)
        gGapXY->Draw("P same");
    c4->cd(2);
    hXZ->DrawClone("colz");
    if(gGapXZ)
        gGapXZ->Draw("P same");
    c4->cd(3);
    hYZ->DrawClone("colz");
    if(gGapYZ)
        gGapYZ->Draw("P same");
}