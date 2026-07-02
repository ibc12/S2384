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

struct ProjectedPoint
{
    double s;
    ROOT::Math::XYZPointF pos3D;
};

// Igual que antes: proyecta los voxeles sobre la línea del cluster, ordenados por s.
std::vector<ProjectedPoint>
GetProjectedPositions(ActRoot::Cluster lightCl, const ROOT::Math::XYZPointF& ref, double driftFactor = 1.)
{
    auto voxels = lightCl.GetRefToVoxels();
    lightCl.ScaleVoxels(2, driftFactor);
    lightCl.ReFit();
    auto line = lightCl.GetRefToLine();
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

// Resumen de gaps de un evento: gap máximo (absoluto y normalizado),
// número de voxeles, y rango total en s.
// La normalización es importante: un gap de 5mm puede ser enorme en una
// traza corta y densa, o insignificante en una traza larga y dispersa.
// Por eso dividimos por el gap "típico" (mediana de los gaps del propio
// evento) para tener una métrica adimensional comparable entre eventos.
struct GapSummary
{
    double maxGap = -1;
    double maxGapNorm = -1; // maxGap / mediana(gaps)
    double medianGap = -1;
    double sRange = -1;
    int nVoxels = 0;
    bool valid = false;
};

GapSummary ComputeGapSummary(const std::vector<ProjectedPoint>& points)
{
    GapSummary res;
    res.nVoxels = (int)points.size();
    if(points.size() < 3) // con <3 puntos no hay "gap típico" fiable
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

    if(res.medianGap > 0)
        res.maxGapNorm = res.maxGap / res.medianGap;

    res.sRange = points.back().s - points.front().s;
    res.valid = true;
    return res;
}

void TracksMissingPadsZProjection()
{
    ROOT::EnableImplicitMT(); // ahora sí, al recorrer todos los eventos conviene multithreadear

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

    // Sin filtro de evento: recorremos todo el chain.
    // Si quieres mantener algún filtro físico (p.ej. solo L1), déjalo aquí.
    auto defSummary = df.Define("gapSummary",
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
                          .Define("nVoxelsG", [](const GapSummary& g) { return g.nVoxels; }, {"gapSummary"})
                          .Define("sRange", [](const GapSummary& g) { return g.sRange; }, {"gapSummary"})
                          .Filter([](const GapSummary& g) { return g.valid; }, {"gapSummary"}); // descarta eventos con <3 voxeles

    // --- Distribución del gap máximo absoluto (mm) ---
    auto hMaxGap = defSummary.Histo1D({"hMaxGap", "Max gap per event;max gap [mm];Events", 200, 0, 50}, "maxGap");

    // --- Distribución del gap máximo normalizado (adimensional) ---
    // Esta es probablemente la más útil para fijar un umbral: un maxGapNorm ~1-2
    // es "espaciado normal", mientras que valores >>1 indican un salto sospechoso
    // independientemente de la longitud/orientación/densidad de la traza.
    auto hMaxGapNorm =
        defSummary.Histo1D({"hMaxGapNorm", "Max gap / median gap;maxGap / medianGap;Events", 200, 0, 50}, "maxGapNorm");

    // --- Correlación: ¿el maxGap depende del número de voxeles (trazas cortas vs largas)? ---
    // Importante para saber si necesitas normalizar por longitud/densidad
    // antes de fijar un umbral global, o si el umbral debe depender de nVoxels.
    auto hGapVsNVoxels = defSummary.Histo2D({"hGapVsNVoxels", "Max gap vs nVoxels;nVoxels;max gap [mm]", 100, 0, 500,
                                             100, 0, 50},
                                            "nVoxelsG", "maxGap");

    auto hGapNormVsNVoxels = defSummary.Histo2D({"hGapNormVsNVoxels", "Max gap norm vs nVoxels;nVoxels;maxGap/medianGap",
                                                 100, 0, 500, 100, 0, 50},
                                                "nVoxelsG", "maxGapNorm");

    // --- Correlación con el rango total de la traza (sRange) ---
    // Traza larga = más oportunidades de tener un gap grande por azar.
    auto hGapVsSRange = defSummary.Histo2D({"hGapVsSRange", "Max gap vs track length;s range [mm];max gap [mm]", 100,
                                            0, 300, 100, 0, 50},
                                           "sRange", "maxGap");

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

    std::cout << "Total events processed: " << *defSummary.Count() << std::endl;
}