// ============================================================================
// PlotTPCEvent_CompareExpSim.C
//
// Para una lista de eventos (run,entry) y una lista de chargeThresholds,
// compara la distribucion de carga por pad, normalizada a su propio maximo
// (Q_pad/Q_max), entre el experimento y la simulacion correspondiente.
//
// Para cada evento y cada threshold se generan:
//   1) Perfil de carga normalizado a lo largo de la traza, exp vs sim
//      superpuestos + panel de ratio sim/exp.
//   2) Proyecciones 2D (XY,XZ,YZ) normalizadas, exp arriba / sim abajo, con
//      la MISMA escala de color [0,1].
// y, por evento, un resumen de RMS(perfil_sim - perfil_exp) vs threshold, mas
// un resumen combinado de esas curvas para los 6 eventos a la vez.
//
// El dato experimental se lee UNA sola vez con RDataFrame, filtrando por
// (fRun,fEntry) para quedarnos solo con los 6 eventos de interes, y
// quedandonos, para cada uno, solo con los clusters fLightIdx/fHeavyIdx de
// ActRoot::MergerData (ver LoadExperimentalVoxelMaps).
//
// !! COSAS A REVISAR ANTES DE CORRERLA !!
//  - expFile / expTree en la seccion de configuracion del main.
//  - Los nombres de cabecera ActMergerData.h / ActTPCData.h si en tu build
//    se llaman distinto.
//  - Se asume fLightIdx/fHeavyIdx == -1 cuando esa particula no se
//    identifico en el evento.
//  - chargeThresholds esta en las mismas unidades (electrones) que la
//    ganancia Gmean de la simulacion, para que el mismo corte tenga sentido
//    aplicado a ambos mapas.
//  - No se activa ROOT::EnableImplicitMT(): el Foreach de RDataFrame corre
//    en un solo hilo, que es lo seguro dado que escribe en un std::map
//    compartido sin proteccion (mutex). Si activas MT, hay que anadir un
//    mutex alrededor de la insercion en 'result'.
// ============================================================================

#include "ActCluster.h"
#include "ActLine.h"
#include "ActMergerData.h"
#include "ActSRIM.h"
#include "ActTPCData.h"
#include "ActTPCParameters.h"
#include "ActVoxel.h"

#include <ROOT/RDataFrame.hxx>
#include <random>

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TRandom3.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <Math/Point3D.h>
#include <Math/Vector3D.h>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "../../PrettyStyle.C"

using XYZPoint = ROOT::Math::XYZPoint;
using XYZVector = ROOT::Math::XYZVector;

// ============================================================
// Geometry
// ============================================================
constexpr float voxelSize = 2.0;            // mm
ActRoot::TPCParameters tpc {"Actar"};       // TPC parameters
constexpr double Gmean = 3000.0;            // Mean gain
constexpr double theta = 0.7;               // Polya parameter
using voxelKey = std::tuple<int, int, int>; // ix,iy,iz

constexpr int yMinExclusionZone = 56;
constexpr int yMaxExclusionZone = 71;
constexpr bool useExclusionZone = true; // poner a false para comparar sin filtrar

// ============================================================
// Polya function for gain distribution in Micromegas
// ============================================================
double Polya(double Gmean, double theta)
{
    static std::mt19937 gen(std::random_device {}());
    std::gamma_distribution<double> gamma(theta + 1.0, Gmean / (theta + 1.0));
    return gamma(gen);
}

// ============================================================
// Add charge to voxel map
// ============================================================
void AddChargeToVoxel(double x, double y, double z, double charge, std::map<voxelKey, ActRoot::Voxel>& voxelMap)
{
    int ix = std::floor(x / voxelSize);
    int iy = std::floor(y / voxelSize);
    int iz = std::floor(z / voxelSize);

    voxelKey key {ix, iy, iz};

    if(ix < 0 || ix >= (tpc.X() / voxelSize))
        return;
    if(iy < 0 || iy >= (tpc.Y() / voxelSize))
        return;

    auto it = voxelMap.find(key);
    if(it == voxelMap.end())
    {
        ActRoot::Voxel v;
        v.SetPosition(ActRoot::Voxel::XYZPointF {float(ix), float(iy), float(iz)});
        v.SetCharge(charge);
        voxelMap.emplace(key, v);
    }
    else
    {
        it->second.SetCharge(it->second.GetCharge() + charge);
    }
}

// ============================================================
// Divide segment -> electrons -> voxels
// ============================================================
void DivideSegmentInPortions(double eLoss, int nPortions, const XYZPoint& center,
                             std::map<voxelKey, ActRoot::Voxel>& voxelMap, std::vector<XYZPoint>& electrons,
                             ActRoot::TPCParameters& tpc)
{
    if(eLoss <= 0 || nPortions <= 0)
        return;

    const double W = 36.4 * 0.95 + 34.3 * 0.05; // eV / electron in 95% D2 + 5% CF4
    const double diffT = 0.06;                  // mm / sqrt(mm)
    const double diffL = 0;                     // mm / sqrt(mm)

    double h = center.Z();
    if(h <= 0)
        return;

    double sigmaT = diffT * std::sqrt(h);
    double sigmaL = diffL * std::sqrt(h);

    const double difXmax = tpc.X() + 20;
    const double difYmax = tpc.Y() + 20;

    if(center.X() < -20 || center.X() > difXmax || center.Y() < -20 || center.Y() > difYmax)
        return;

    double portionE = eLoss / nPortions; // MeV
    double meanNe = portionE * 1e6 / W;

    for(int i = 0; i < nPortions; i++)
    {
        int nElectrons = gRandom->Poisson(meanNe);
        if(nElectrons <= 0)
            continue;

        for(int e = 0; e < nElectrons; e++)
        {
            double x = gRandom->Gaus(center.X(), sigmaT);
            double y = gRandom->Gaus(center.Y(), sigmaT);
            double z = gRandom->Gaus(center.Z(), sigmaL);

            double gain = Polya(Gmean, theta);
            AddChargeToVoxel(x, y, z, gain, voxelMap);
            electrons.emplace_back(x, y, z);
        }
    }
}

// ============================================================
// Divide track using SRIM
// ============================================================
void DivideTrackInSegments(ActPhysics::SRIM* srim, double range, const XYZVector& dirIn, const XYZPoint& rp,
                           double step, int nSub, std::map<voxelKey, ActRoot::Voxel>& voxelMap,
                           std::vector<XYZPoint>& electrons, ActRoot::TPCParameters& tpc,
                           std::string particleType = "heavy")
{
    XYZVector dir = dirIn.Unit();
    double E = srim->EvalInverse(particleType, range);

    for(double r = 0; r < range; r += step)
    {
        double Epost = srim->Slow(particleType, E, step);
        double eLoss = E - Epost;
        E = Epost;

        if(eLoss <= 0)
            continue;

        XYZPoint center = rp + dir * (r + 0.5 * step);
        if(center.Z() <= 0 || center.Z() > 256)
            continue;

        DivideSegmentInPortions(eLoss, nSub, center, voxelMap, electrons, tpc);
    }
}

// ============================================================
// Build profile (TGraph subdividido + TH1D) a partir de un mapa de voxeles.
// ============================================================
std::pair<TGraph*, TH1D*> GetChargeProfile(const std::map<voxelKey, ActRoot::Voxel>& voxelMap, bool subdivideVoxels)
{
    std::vector<ActRoot::Voxel> voxels;
    voxels.reserve(voxelMap.size());
    for(const auto& [k, v] : voxelMap)
        voxels.push_back(v);

    ActRoot::Cluster cluster;
    cluster.SetVoxels(voxels);
    cluster.ReFit();

    ActRoot::Line line = cluster.GetLine();
    line.Scale(voxelSize, voxelSize);
    auto p0 = line.GetPoint();
    auto u = line.GetDirection().Unit();

    double sMin = 1e9;
    double sMax = -1e9;

    for(auto v : voxels)
    {
        auto pos = v.GetPosition();
        v.SetPosition({pos.X() * voxelSize, pos.Y() * voxelSize, pos.Z() * voxelSize});
        XYZVector d = XYZVector(v.GetPosition()) - XYZVector(p0);
        double s = d.Dot(u);
        sMin = std::min(sMin, s);
        sMax = std::max(sMax, s);
    }

    double length = sMax - sMin;
    if(length <= 0)
        length = voxelSize;

    double ds = voxelSize * std::max({std::abs(u.X()), std::abs(u.Y()), std::abs(u.Z())});
    ds = std::max(ds, 0.5 * voxelSize);
    int nBins = std::max(1, int(std::ceil((length + 10) / ds)));

    auto* graph = new TGraph();
    graph->SetName("chargeProfileG");
    graph->SetTitle("Charge profile;Track length [mm];Q");

    auto* hist = new TH1D("chargeProfileH", "Charge profile;Track length [mm];Q", nBins, sMin - 5, sMax + 5);

    int nDiv = 3;
    double div = 1.0 / nDiv;

    for(auto v : voxels)
    {
        auto pos = v.GetPosition();
        double q = v.GetCharge();

        if(!subdivideVoxels)
        {
            v.SetPosition({(pos.X() + 0.5f) * voxelSize, (pos.Y() + 0.5f) * voxelSize, (pos.Z() + 0.5f) * voxelSize});
            XYZVector d = XYZVector(v.GetPosition()) - XYZVector(p0);
            double s = d.Dot(u);
            graph->SetPoint(graph->GetN(), s, q);
            hist->Fill(s, q);
        }
        else
        {
            double qSub = q / (nDiv * nDiv * nDiv);
            for(int ix = -1; ix <= 1; ix++)
                for(int iy = -1; iy <= 1; iy++)
                    for(int iz = -1; iz <= 1; iz++)
                    {
                        XYZPoint miniPos(pos.X() * voxelSize + ix * div * voxelSize,
                                         pos.Y() * voxelSize + iy * div * voxelSize,
                                         pos.Z() * voxelSize + iz * div * voxelSize);
                        XYZVector d = XYZVector(miniPos) - XYZVector(p0);
                        double s = d.Dot(u);
                        graph->SetPoint(graph->GetN(), s, qSub);
                        hist->Fill(s, qSub);
                    }
        }
    }

    return {graph, hist};
}

// ============================================================
// Suma dos mapas de voxeles, acumulando la carga en los pads que coinciden.
// ============================================================
std::map<voxelKey, ActRoot::Voxel>
MergeVoxelMaps(const std::map<voxelKey, ActRoot::Voxel>& mapA, const std::map<voxelKey, ActRoot::Voxel>& mapB)
{
    std::map<voxelKey, ActRoot::Voxel> merged = mapA;
    for(const auto& [key, v] : mapB)
    {
        auto it = merged.find(key);
        if(it == merged.end())
            merged.emplace(key, v);
        else
            it->second.SetCharge(it->second.GetCharge() + v.GetCharge());
    }
    return merged;
}

// ============================================================
// Filtra los pads dentro de la exclusion zone en iy.
// ============================================================
std::map<voxelKey, ActRoot::Voxel> ApplyExclusionZone(const std::map<voxelKey, ActRoot::Voxel>& voxelMap)
{
    if(!useExclusionZone)
        return voxelMap;

    std::map<voxelKey, ActRoot::Voxel> filtered;
    for(const auto& [key, v] : voxelMap)
    {
        int iy = std::get<1>(key);
        if(iy < yMinExclusionZone || iy > yMaxExclusionZone)
            filtered.emplace(key, v);
    }
    return filtered;
}

// ============================================================
// Se queda solo con los pads cuya carga (cruda, en electrones) supera el
// threshold dado. Se aplica IGUAL a exp y sim antes de normalizar.
// ============================================================
std::map<voxelKey, ActRoot::Voxel>
FilterByThreshold(const std::map<voxelKey, ActRoot::Voxel>& voxelMap, double threshold)
{
    std::map<voxelKey, ActRoot::Voxel> filtered;
    for(const auto& [key, v] : voxelMap)
        if(v.GetCharge() > threshold)
            filtered.emplace(key, v);
    return filtered;
}

// ============================================================
// Maximo de carga de un mapa (0 si esta vacio)
// ============================================================
double GetMaxCharge(const std::map<voxelKey, ActRoot::Voxel>& voxelMap)
{
    double qMax = 0;
    for(const auto& [key, v] : voxelMap)
    {
        if(v.GetCharge() > qMax)
            qMax = v.GetCharge();
    }
    return qMax;
}

// ============================================================
// Devuelve una COPIA del mapa con la carga dividida por qMax (Q -> [0,1])
// ============================================================
std::map<voxelKey, ActRoot::Voxel> NormalizeVoxelMap(const std::map<voxelKey, ActRoot::Voxel>& voxelMap, double qMax)
{
    std::map<voxelKey, ActRoot::Voxel> norm = voxelMap;
    if(qMax <= 0)
        return norm;
    for(auto& [key, v] : norm)
        v.SetCharge(v.GetCharge() / qMax);
    return norm;
}

// ============================================================
// Proyecciones 2D (XY, XZ, YZ) de un mapa YA normalizado, escala fija [0,1].
// ============================================================
struct ProjectionPlots
{
    TH2D* hXY = nullptr;
    TH2D* hXZ = nullptr;
    TH2D* hYZ = nullptr;
};

ProjectionPlots BuildProjectionPlots(const std::map<voxelKey, ActRoot::Voxel>& normMap, ActRoot::TPCParameters& tpc,
                                     const TString& tag, const TString& label)
{
    ProjectionPlots p;

    p.hXY = new TH2D("hXYq_" + tag, label + " - XY;X [pad];Y [pad]", tpc.X() / voxelSize, 0, tpc.X() / voxelSize,
                     tpc.Y() / voxelSize, 0, tpc.Y() / voxelSize);
    p.hXZ = new TH2D("hXZq_" + tag, label + " - XZ;X [pad];Z [pad]", tpc.X() / voxelSize, 0, tpc.X() / voxelSize,
                     tpc.Z() / voxelSize, 0, tpc.Z() / voxelSize);
    p.hYZ = new TH2D("hYZq_" + tag, label + " - YZ;Z [pad];Y [pad]", tpc.Z() / voxelSize, 0, tpc.Z() / voxelSize,
                     tpc.Y() / voxelSize, 0, tpc.Y() / voxelSize);

    for(const auto& [key, v] : normMap)
    {
        const auto& pos = v.GetPosition();
        double q = v.GetCharge();
        p.hXY->Fill(pos.X(), pos.Y(), q);
        p.hXZ->Fill(pos.X(), pos.Z(), q);
        p.hYZ->Fill(pos.Z(), pos.Y(), q);
    }

    for(auto* h : {p.hXY, p.hXZ, p.hYZ})
    {
        h->SetMinimum(0);
        h->SetMaximum(1);
    }

    return p;
}

// ============================================================
// Rellena un histograma NUEVO, con un binning comun, a partir de los puntos
// (x,y) de un TGraph, aplicando un desplazamiento en x.
// ============================================================
TH1D* FillCommonProfile(TGraph* g, double shift, int nBins, double xmin, double xmax, const TString& name,
                        const TString& title)
{
    auto* h = new TH1D(name, title, nBins, xmin, xmax);
    double* xs = g->GetX();
    double* ys = g->GetY();
    for(int i = 0; i < g->GetN(); i++)
        h->Fill(xs[i] + shift, ys[i]);
    return h;
}

// ============================================================
// Lee UNA sola vez, con RDataFrame, los mapas de voxeles experimentales de
// TODOS los eventos (run,entry) pedidos, quedandose solo con los clusters
// fLightIdx/fHeavyIdx de MergerData (igual que simulas solo light+heavy).
// Devuelve un mapa {(run,entry) -> voxelMap}; si un evento no aparece en el
// arbol, simplemente no tendra entrada (comprobar con .find en el caller).
// ============================================================
std::map<std::pair<int, int>, std::map<voxelKey, ActRoot::Voxel>>
LoadExperimentalVoxelMaps(const std::string& file, const std::string& treeName, const std::vector<int>& runs,
                          const std::vector<int>& entries)
{
    std::map<std::pair<int, int>, std::map<voxelKey, ActRoot::Voxel>> result;

    std::set<std::pair<int, int>> wanted;
    for(size_t i = 0; i < runs.size(); i++)
        wanted.insert({runs[i], entries[i]});

    ROOT::RDataFrame df {treeName.c_str(), file.c_str()};

    auto dfFiltered =
        df.Filter([wanted](ActRoot::MergerData& m) { return wanted.count({m.fRun, m.fEntry}) > 0; }, {"MergerData"});

    dfFiltered.Foreach(
        [&result](ActRoot::MergerData& m, ActRoot::TPCData& tpcData)
        {
            std::map<voxelKey, ActRoot::Voxel> voxelMap;

            std::vector<int> idxs;
            if(m.fLightIdx >= 0)
                idxs.push_back(m.fLightIdx);
            if(m.fHeavyIdx >= 0)
                idxs.push_back(m.fHeavyIdx);

            for(int idx : idxs)
            {
                if(idx < 0 || idx >= (int)tpcData.fClusters.size())
                    continue;

                auto& cluster = tpcData.fClusters[idx];
                for(auto& voxel : cluster.GetRefToVoxels())
                {
                    auto pos = voxel.GetPosition();
                    voxelKey key {(int)pos.X(), (int)pos.Y(), (int)pos.Z()};
                    double q = voxel.GetCharge();

                    auto it = voxelMap.find(key);
                    if(it == voxelMap.end())
                        voxelMap.emplace(key, voxel);
                    else
                        it->second.SetCharge(it->second.GetCharge() + q);
                }
            }

            result[{m.fRun, m.fEntry}] = std::move(voxelMap);
        },
        {"MergerData", "TPCData"});

    std::cout << "[LoadExperimentalVoxelMaps] Encontrados " << result.size() << " / " << runs.size()
              << " eventos pedidos en " << file << std::endl;

    return result;
}

// ============================================================
// MAIN
// ============================================================
void PlotTPCEvent_ChargeDistribution()
{
    PrettyStyle();
    gRandom->SetSeed(0);
    TH1::AddDirectory(kFALSE);

    // ================= 1. SRIM tables =================
    auto* srim = new ActPhysics::SRIM;
    srim->ReadTable("light", "../../Calibrations/SRIM/1H_900mb_CF4_95-5.txt");
    srim->ReadTable("lightD", "../../Calibrations/SRIM/2H_900mb_CF4_95-5.txt");
    srim->ReadTable("lightT", "../../Calibrations/SRIM/3H_900mb_CF4_95-5.txt");
    srim->ReadTable("heavy", "../../Calibrations/SRIM/11Li_900mb_CF4_95-5.txt");

    // ================= 2. Events =================
    // std::vector<std::string> eventNames = {"r66e292", "r66e973", "r66e1554", "r66e6583", "r66e8012", "r67e23132"};
    // std::vector<int> runs = {66, 66, 66, 66, 66, 67};
    // std::vector<int> entries = {292, 973, 1554, 6583, 8012, 23132};
    //
    // std::vector<double> thetaLights = {79, 77.2, 81.14, 149.5, 128.5, 149.4};
    // std::vector<double> phiLights = {-122.8, 105.5, -100.6, -148.3, -87, -44};
    // std::vector<double> thetaHeavy = {4.8, 5.5, 3.64, 2.8, 2.76, 1.36};
    // std::vector<double> phiHeavy = {43.43, -72.73, 72.3, 30.57, 104.3, 96.3};
    // std::vector<double> rp_x = {200.5, 16.4, 198., 208., 129.7, 149.8};
    // std::vector<double> rp_y = {128.5, 123.6, 126., 126., 123.6, 123.4};
    // std::vector<double> ranges = {30.8, 115.8, 23.4, 99.23, 136, 112.7};
    // std::vector<std::string> lightStrings = {"lightD", "lightD", "lightD", "light", "light", "light"};

    // ================= 2. Events without short ones =================
    std::vector<std::string> eventNames = {"r66e152", "r66e973", "r66e1144", "r66e6583", "r66e8012", "r67e23132"};
    std::vector<int> runs = {66, 66, 66, 66, 66, 67};
    std::vector<int> entries = {152, 973, 1144, 6583, 8012, 23132};

    std::vector<double> thetaLights = {79.43, 77.2, 78.3, 149.5, 128.5, 149.4};
    std::vector<double> phiLights = {58.5, 105.5, -101.7, -148.3, -87, -44};
    std::vector<double> thetaHeavy = {4.44, 5.5, 5.73, 2.8, 2.76, 1.36};
    std::vector<double> phiHeavy = {-118.5, -72.73, 66.6, 30.57, 104.3, 96.3};
    std::vector<double> rp_x = {138.9, 16.4, 202.4, 208., 129.7, 149.8};
    std::vector<double> rp_y = {125., 123.6, 124.1, 126., 123.6, 123.4};
    std::vector<double> ranges = {59.19, 115.8, 82.4, 99.23, 136, 112.7};
    std::vector<std::string> lightStrings = {"lightD", "lightD", "lightD", "light", "light", "light"};

    // IMPORTANT: use a threshold appropriate for the experimental charge scale.
    // Charge thresholds to compare. Adjust these values to the experimental charge scale.
    std::vector<double> chargeThresholds = {1e2, 7e2, 1e3, 5e3, 7e3, 9e3, 1e4, 2e4, 3e4, 4e4, 7e4, 1e5, 4e5, 1e6};

    // ================= 3. Experimental data =================
    std::string expFile = "../../PostAnalysis/Outputs/tree_preprocess_F_7Li.root";
    std::string expTree = "PreProcessed_Tree";
    auto expMapsByEvent = LoadExperimentalVoxelMaps(expFile, expTree, runs, entries);

    const size_t nThr = chargeThresholds.size();

    // Results grouped by threshold. The event index is preserved in every vector.
    std::vector<std::vector<TH1D*>> profilesSimByThr(nThr);
    std::vector<std::vector<TH1D*>> profilesExpByThr(nThr);
    std::vector<std::vector<TH1D*>> ratiosByThr(nThr);
    std::vector<std::vector<double>> rmsByThr(nThr);
    std::vector<std::vector<std::string>> eventNamesByThr(nThr);

    // ================= Loop over events =================
    for(size_t i = 0; i < eventNames.size(); i++)
    {
        const std::string& eventName = eventNames[i];
        std::cout << "\n\n########## " << eventName << " (run " << runs[i] << ", entry " << entries[i] << ") ##########"
                  << std::endl;

        // ---- Simulation of this event ----
        double thLight = thetaLights[i] * TMath::DegToRad();
        double phLight = phiLights[i] * TMath::DegToRad();
        XYZVector dirLight(std::cos(thLight), std::sin(thLight) * std::sin(phLight),
                           std::sin(thLight) * std::cos(phLight));

        double thHeavy = thetaHeavy[i] * TMath::DegToRad();
        double phHeavy = phiHeavy[i] * TMath::DegToRad();
        XYZVector dirHeavy(std::cos(thHeavy), std::sin(thHeavy) * std::sin(phHeavy),
                           std::sin(thHeavy) * std::cos(phHeavy));

        XYZPoint rp(rp_x[i], rp_y[i], 110);

        std::map<voxelKey, ActRoot::Voxel> voxelMapLight, voxelMapHeavy;
        std::vector<XYZPoint> electronsLight, electronsHeavy;

        DivideTrackInSegments(srim, ranges[i], dirLight, rp, 2, 5, voxelMapLight, electronsLight, tpc, lightStrings[i]);
        DivideTrackInSegments(srim, 3000, dirHeavy, rp, 2, 5, voxelMapHeavy, electronsHeavy, tpc);

        auto voxelMapSimRaw = ApplyExclusionZone(MergeVoxelMaps(voxelMapLight, voxelMapHeavy));

        // ---- Experimental event ----
        std::map<voxelKey, ActRoot::Voxel> voxelMapExpRaw;
        bool hasExp = false;

        auto itExp = expMapsByEvent.find({runs[i], entries[i]});
        if(itExp != expMapsByEvent.end() && !itExp->second.empty())
        {
            voxelMapExpRaw = ApplyExclusionZone(itExp->second);
            hasExp = !voxelMapExpRaw.empty();
        }

        if(!hasExp)
            std::cerr << "[AVISO] Sin dato experimental utilizable para " << eventName << std::endl;

        // ================= Loop over thresholds =================
        for(size_t k = 0; k < nThr; k++)
        {
            const double thr = chargeThresholds[k];

            auto simThr = FilterByThreshold(voxelMapSimRaw, thr);
            if(simThr.empty())
            {
                std::cout << "  [thr=" << thr << "] mapa simulado vacio tras el corte." << std::endl;
                continue;
            }

            bool hasExpThr = false;
            std::map<voxelKey, ActRoot::Voxel> expThr;
            if(hasExp)
            {
                expThr = FilterByThreshold(voxelMapExpRaw, 0); // Exp already has threshold applied in the preprocessing
                hasExpThr = !expThr.empty();
            }

            const double qMaxSim = GetMaxCharge(simThr);
            auto normSim = NormalizeVoxelMap(simThr, qMaxSim);

            double qMaxExp = 0;
            std::map<voxelKey, ActRoot::Voxel> normExp;
            if(hasExpThr)
            {
                qMaxExp = GetMaxCharge(expThr);
                normExp = NormalizeVoxelMap(expThr, qMaxExp);
            }

            std::cout << "  [thr=" << thr << " e-] Qmax_sim=" << qMaxSim
                      << (hasExpThr ? Form(" Qmax_exp=%g", qMaxExp) : " (sin exp tras el corte)") << std::endl;

            // ---- Charge profiles ----
            auto [graphSim, histSim] = GetChargeProfile(normSim, true);
            const double shiftSim = -histSim->GetXaxis()->GetXmin();
            const double lengthSim = histSim->GetXaxis()->GetXmax() - histSim->GetXaxis()->GetXmin();

            TGraph* graphExpP = nullptr;
            TH1D* histExpP = nullptr;
            double shiftExp = 0;
            double lengthExp = 0;

            if(hasExpThr)
            {
                std::tie(graphExpP, histExpP) = GetChargeProfile(normExp, true);
                shiftExp = -histExpP->GetXaxis()->GetXmin();
                lengthExp = histExpP->GetXaxis()->GetXmax() - histExpP->GetXaxis()->GetXmin();
            }

            const double sMaxCommon = std::max(lengthSim, lengthExp) + 5;
            const int nBinsCommon = std::max(1, static_cast<int>(std::ceil(sMaxCommon / voxelSize)));

            TH1D* hSimC = FillCommonProfile(graphSim, shiftSim, nBinsCommon, 0, sMaxCommon,
                                            TString::Format("hProfileSim_%s_thr%zu", eventName.c_str(), k),
                                            TString::Format("Normalized charge profile - %s (thr=%.3g e^{-});"
                                                            "Track length from start [mm];Q/Q_{max}",
                                                            eventName.c_str(), thr));

            hSimC->SetDirectory(nullptr);
            hSimC->SetLineColor(kRed + 1);
            hSimC->SetLineWidth(2);
            hSimC->SetFillColorAlpha(kRed + 1, 0.25);

            TH1D* hExpC = nullptr;
            if(hasExpThr)
            {
                hExpC = FillCommonProfile(graphExpP, shiftExp, nBinsCommon, 0, sMaxCommon,
                                          TString::Format("hProfileExp_%s_thr%zu", eventName.c_str(), k),
                                          hSimC->GetTitle());

                hExpC->SetDirectory(nullptr);
                hExpC->SetLineColor(kBlack);
                hExpC->SetMarkerColor(kBlack);
                hExpC->SetMarkerStyle(20);
                hExpC->SetMarkerSize(0.7);
                hExpC->SetLineWidth(2);
            }

            // Store profiles. Keep nullptr for experiment when unavailable.
            profilesSimByThr[k].push_back(hSimC);
            profilesExpByThr[k].push_back(hExpC);
            eventNamesByThr[k].push_back(eventName);

            // ---- RMS and ratio ----
            double rms = -1;
            TH1D* hRatio = nullptr;

            if(hasExpThr && hExpC != nullptr)
            {
                double sumSq = 0;
                int nComp = 0;

                for(int b = 1; b <= hSimC->GetNbinsX(); b++)
                {
                    const double qs = hSimC->GetBinContent(b);
                    const double qe = hExpC->GetBinContent(b);

                    if(qs == 0 && qe == 0)
                        continue;

                    sumSq += (qs - qe) * (qs - qe);
                    nComp++;
                }

                if(nComp > 0)
                    rms = std::sqrt(sumSq / nComp);

                hRatio = static_cast<TH1D*>(hSimC->Clone(TString::Format("hRatio_%s_thr%zu", eventName.c_str(), k)));
                hRatio->SetDirectory(nullptr);

                // Manual division avoids NaN/Inf bins when Exp = 0.
                for(int b = 1; b <= hRatio->GetNbinsX(); b++)
                {
                    const double qs = hSimC->GetBinContent(b);
                    const double qe = hExpC->GetBinContent(b);

                    if(qe != 0)
                        hRatio->SetBinContent(b, qs / qe);
                    else
                        hRatio->SetBinContent(b, 0);

                    hRatio->SetBinError(b, 0);
                }

                hRatio->SetTitle(TString::Format("Sim/Exp ratio - %s (thr=%.3g e^{-});"
                                                 "Track length from start [mm];Q_{sim}/Q_{exp}",
                                                 eventName.c_str(), thr));
                hRatio->SetLineColor(kAzure + 2);
                hRatio->SetMarkerColor(kAzure + 2);
                hRatio->SetMarkerStyle(20);
                hRatio->SetMarkerSize(0.7);
                hRatio->SetLineWidth(2);
            }

            rmsByThr[k].push_back(rms);
            ratiosByThr[k].push_back(hRatio);

            // ---- 2D normalized projections: keep individual canvas ----
            auto projSim = BuildProjectionPlots(normSim, tpc, TString::Format("sim_%s_thr%zu", eventName.c_str(), k),
                                                "Simulation");

            TCanvas* cProj = new TCanvas(TString::Format("cProj_%s_thr%zu", eventName.c_str(), k),
                                         TString::Format("Normalized maps %s thr=%.3g", eventName.c_str(), thr), 1500,
                                         hasExpThr ? 900 : 500);

            if(hasExpThr)
            {
                auto projExp = BuildProjectionPlots(
                    normExp, tpc, TString::Format("exp_%s_thr%zu", eventName.c_str(), k), "Experiment");

                cProj->Divide(3, 2);
                cProj->cd(1);
                projExp.hXY->Draw("COLZ");
                cProj->cd(2);
                projExp.hXZ->Draw("COLZ");
                cProj->cd(3);
                projExp.hYZ->Draw("COLZ");
                cProj->cd(4);
                projSim.hXY->Draw("COLZ");
                cProj->cd(5);
                projSim.hXZ->Draw("COLZ");
                cProj->cd(6);
                projSim.hYZ->Draw("COLZ");
            }
            else
            {
                cProj->Divide(3, 1);
                cProj->cd(1);
                projSim.hXY->Draw("COLZ");
                cProj->cd(2);
                projSim.hXZ->Draw("COLZ");
                cProj->cd(3);
                projSim.hYZ->Draw("COLZ");
            }
        }
    }

    // ============================================================
    // COMBINED CANVASES: one set for each threshold
    // ============================================================
    for(size_t k = 0; k < nThr; k++)
    {
        const double thr = chargeThresholds[k];
        const int nEvents = static_cast<int>(eventNamesByThr[k].size());
        if(nEvents == 0)
            continue;

        const int nCols = 3;
        const int nRows = std::max(1, static_cast<int>(std::ceil(nEvents / 3.0)));

        // ========================================================
        // 1. ALL CHARGE PROFILES
        // ========================================================
        TCanvas* cProfilesAll =
            new TCanvas(TString::Format("cProfilesAll_thr%zu", k),
                        TString::Format("All charge profiles - threshold %.3g e-", thr), 1800, 1000);

        cProfilesAll->Divide(nCols, nRows);

        for(int ievt = 0; ievt < nEvents; ievt++)
        {
            cProfilesAll->cd(ievt + 1);

            TH1D* hSim = profilesSimByThr[k][ievt];
            TH1D* hExp = profilesExpByThr[k][ievt];
            if(!hSim)
                continue;

            hSim->SetTitle(
                TString::Format("%s;Track length from start [mm];Q/Q_{max}", eventNamesByThr[k][ievt].c_str()));
            hSim->Draw("HIST");

            if(hExp)
                hExp->Draw("PE SAME");

            auto* leg = new TLegend(0.58, 0.72, 0.89, 0.89);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(hSim, "Simulation", "f");
            if(hExp)
                leg->AddEntry(hExp, "Experiment", "lp");
            leg->Draw();
        }

        // ========================================================
        // 2. RMS FOR ALL EVENTS AT THIS THRESHOLD
        // ========================================================
        TCanvas* cRmsAll = new TCanvas(TString::Format("cRmsAll_thr%zu", k),
                                       TString::Format("RMS - threshold %.3g e-", thr), 1100, 650);

        auto* hRms = new TH1D(TString::Format("hRmsAll_thr%zu", k),
                              TString::Format("RMS comparison - threshold %.3g e-;Event;"
                                              "RMS(Q/Q_{max,sim} - Q/Q_{max,exp})",
                                              thr),
                              nEvents, 0, nEvents);

        double maxRms = 0;
        bool hasAnyRms = false;

        for(int ievt = 0; ievt < nEvents; ievt++)
        {
            const double rms = rmsByThr[k][ievt];
            hRms->GetXaxis()->SetBinLabel(ievt + 1, eventNamesByThr[k][ievt].c_str());

            if(rms >= 0)
            {
                hRms->SetBinContent(ievt + 1, rms);
                maxRms = std::max(maxRms, rms);
                hasAnyRms = true;
            }
        }

        hRms->SetLineWidth(2);
        hRms->SetMarkerStyle(20);
        hRms->SetMarkerSize(1.2);
        hRms->SetMinimum(0);
        hRms->SetMaximum(hasAnyRms ? std::max(0.05, 1.2 * maxRms) : 1.0);
        hRms->LabelsOption("v", "X");
        hRms->Draw("HIST P");

        // ========================================================
        // 3. ALL SIM/EXP RATIOS
        // ========================================================
        TCanvas* cRatiosAll = new TCanvas(TString::Format("cRatiosAll_thr%zu", k),
                                          TString::Format("All Sim/Exp ratios - threshold %.3g e-", thr), 1800, 1000);

        cRatiosAll->Divide(nCols, nRows);

        for(int ievt = 0; ievt < nEvents; ievt++)
        {
            cRatiosAll->cd(ievt + 1);

            TH1D* hRatio = ratiosByThr[k][ievt];
            if(!hRatio)
            {
                auto* text = new TPaveText(0.20, 0.40, 0.80, 0.60, "NDC");
                text->SetBorderSize(0);
                text->SetFillStyle(0);
                text->SetTextAlign(22);
                text->AddText(eventNamesByThr[k][ievt].c_str());
                text->AddText("No experimental data after threshold");
                text->Draw();
                continue;
            }

            hRatio->SetTitle(
                TString::Format("%s;Track length from start [mm];Q_{sim}/Q_{exp}", eventNamesByThr[k][ievt].c_str()));

            double ratioMax = 0;
            for(int b = 1; b <= hRatio->GetNbinsX(); b++)
            {
                const double r = hRatio->GetBinContent(b);
                if(r > 0 && std::isfinite(r))
                    ratioMax = std::max(ratioMax, r);
            }

            hRatio->SetMinimum(0);
            hRatio->SetMaximum(ratioMax > 0 ? std::max(2.0, 1.15 * ratioMax) : 2.0);
            hRatio->Draw("PE");

            const double xmin = hRatio->GetXaxis()->GetXmin();
            const double xmax = hRatio->GetXaxis()->GetXmax();
            auto* lineOne = new TLine(xmin, 1, xmax, 1);
            lineOne->SetLineStyle(2);
            lineOne->SetLineColor(kGray + 2);
            lineOne->SetLineWidth(2);
            lineOne->Draw();
        }
    }

    // ============================================================
    // GLOBAL SUMMARY: RMS comparison for all charge thresholds
    // ============================================================
    // One pad per threshold. Each pad shows the RMS obtained for all
    // events at that threshold. The mean RMS over valid events is shown
    // in the pad, making it easy to identify the best threshold.
    {
        const int nSummaryCols = 2;
        const int nSummaryRows = std::max(1, static_cast<int>(std::ceil(nThr / static_cast<double>(nSummaryCols))));

        TCanvas* cRmsSummary = new TCanvas("cRmsSummaryAllThresholds", "RMS comparison for all charge thresholds", 1400,
                                           500 * nSummaryRows);

        cRmsSummary->Divide(nSummaryCols, nSummaryRows);

        std::vector<double> meanRmsByThr(nThr, -1.0);
        std::vector<int> nValidRmsByThr(nThr, 0);

        double globalMaxRms = 0.0;
        for(size_t k = 0; k < nThr; k++)
        {
            for(double rms : rmsByThr[k])
                if(rms >= 0 && std::isfinite(rms))
                    globalMaxRms = std::max(globalMaxRms, rms);
        }

        for(size_t k = 0; k < nThr; k++)
        {
            cRmsSummary->cd(k + 1);
            gPad->SetBottomMargin(0.22);

            const int nEvents = static_cast<int>(eventNamesByThr[k].size());
            if(nEvents == 0)
                continue;

            const double thr = chargeThresholds[k];

            auto* hSummary = new TH1D(TString::Format("hRmsSummary_thr%zu", k),
                                      TString::Format("Threshold = %.3g e^{-};Event;RMS", thr), nEvents, 0, nEvents);

            double sumRms = 0.0;
            int nValid = 0;

            for(int ievt = 0; ievt < nEvents; ievt++)
            {
                hSummary->GetXaxis()->SetBinLabel(ievt + 1, eventNamesByThr[k][ievt].c_str());

                const double rms = rmsByThr[k][ievt];
                if(rms >= 0 && std::isfinite(rms))
                {
                    hSummary->SetBinContent(ievt + 1, rms);
                    sumRms += rms;
                    nValid++;
                }
            }

            if(nValid > 0)
            {
                meanRmsByThr[k] = sumRms / nValid;
                nValidRmsByThr[k] = nValid;
            }

            hSummary->SetMinimum(0.0);
            hSummary->SetMaximum(globalMaxRms > 0 ? 1.15 * globalMaxRms : 1.0);
            hSummary->SetLineWidth(2);
            hSummary->SetMarkerStyle(20);
            hSummary->SetMarkerSize(1.0);
            hSummary->LabelsOption("v", "X");
            hSummary->Draw("HIST P");

            if(nValid > 0)
            {
                auto* meanLine = new TLine(0.0, meanRmsByThr[k], static_cast<double>(nEvents), meanRmsByThr[k]);
                meanLine->SetLineColor(kRed + 1);
                meanLine->SetLineStyle(2);
                meanLine->SetLineWidth(2);
                meanLine->Draw();

                auto* info = new TPaveText(0.58, 0.74, 0.89, 0.89, "NDC");
                info->SetFillStyle(0);
                info->SetBorderSize(1);
                info->SetTextAlign(12);
                info->SetTextSize(0.040);
                info->AddText(TString::Format("Mean RMS = %.4f", meanRmsByThr[k]));
                info->AddText(TString::Format("Valid events = %d/%d", nValid, nEvents));
                info->Draw();
            }
        }

        // ========================================================
        // MEAN RMS VS CHARGE THRESHOLD
        // ========================================================
        // This is the most direct plot for selecting the threshold:
        // the preferred value is the one with the minimum mean RMS,
        // provided it is based on a sufficient number of valid events.
        auto* gMeanRmsVsThr = new TGraph();
        gMeanRmsVsThr->SetName("gMeanRmsVsChargeThreshold");
        gMeanRmsVsThr->SetTitle("Mean RMS vs charge threshold;Charge threshold [e^{-}];Mean RMS");

        double bestThr = -1.0;
        double bestMeanRms = 1e300;
        int bestNValid = 0;

        for(size_t k = 0; k < nThr; k++)
        {
            if(meanRmsByThr[k] < 0)
                continue;

            const int p = gMeanRmsVsThr->GetN();
            gMeanRmsVsThr->SetPoint(p, chargeThresholds[k], meanRmsByThr[k]);

            if(meanRmsByThr[k] < bestMeanRms)
            {
                bestMeanRms = meanRmsByThr[k];
                bestThr = chargeThresholds[k];
                bestNValid = nValidRmsByThr[k];
            }
        }

        if(gMeanRmsVsThr->GetN() > 0)
        {
            TCanvas* cMeanRmsVsThr = new TCanvas("cMeanRmsVsChargeThreshold", "Mean RMS vs charge threshold", 900, 650);

            cMeanRmsVsThr->SetLogx();

            gMeanRmsVsThr->SetMarkerStyle(20);
            gMeanRmsVsThr->SetMarkerSize(1.2);
            gMeanRmsVsThr->SetMarkerColor(kAzure + 2);
            gMeanRmsVsThr->SetLineColor(kAzure + 2);
            gMeanRmsVsThr->SetLineWidth(2);
            gMeanRmsVsThr->Draw("APL");

            if(bestThr > 0)
            {
                auto* bestLine = new TLine(bestThr, cMeanRmsVsThr->GetUymin(), bestThr, cMeanRmsVsThr->GetUymax());
                bestLine->SetLineColor(kRed + 1);
                bestLine->SetLineStyle(2);
                bestLine->SetLineWidth(2);
                bestLine->Draw();

                auto* bestInfo = new TPaveText(0.52, 0.76, 0.89, 0.89, "NDC");
                bestInfo->SetFillStyle(0);
                bestInfo->SetBorderSize(1);
                bestInfo->SetTextAlign(12);
                bestInfo->SetTextSize(0.040);
                bestInfo->AddText("Minimum mean RMS:");
                bestInfo->AddText(TString::Format("Threshold = %.3g e^{-}", bestThr));
                bestInfo->AddText(TString::Format("Mean RMS = %.4f", bestMeanRms));
                bestInfo->AddText(TString::Format("Valid events = %d", bestNValid));
                bestInfo->Draw();

                std::cout << "\n============================================================\n";
                std::cout << "BEST CHARGE THRESHOLD ACCORDING TO MEAN RMS\n";
                std::cout << "  Threshold = " << bestThr << " e-\n";
                std::cout << "  Mean RMS  = " << bestMeanRms << "\n";
                std::cout << "  Valid events = " << bestNValid << "\n";
                std::cout << "============================================================\n";
            }
        }
    }
}
