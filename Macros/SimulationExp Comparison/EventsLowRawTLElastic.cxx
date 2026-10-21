#include "ActMergerData.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "../../PrettyStyle.C"

// ---------------- Configuración ----------------
constexpr double rawTLCut = 40.;  // mm: RawTL < rawTLCut -> muestra "low"
constexpr double padSize = 2.;    // mm per pad
constexpr double detMax = 256.;   // mm
constexpr double thetaMin = 130.; // deg
constexpr double thetaMax = 170.; // deg
constexpr int thetaBins = 40;
constexpr double eVertexMax = 3.; // MeV
constexpr double qMax = 20000.;   // Rango de fQtotal (AJUSTAR a tus unidades)

// Posición (X, Y) en pads del último vóxel del cluster ligero
std::pair<double, double> GetLastVoxelXY(ActRoot::MergerData& m, ActRoot::TPCData& tpc)
{
    auto lightCl = tpc.fClusters[m.fLightIdx]; // copia: la modificamos
    auto line = lightCl.GetLine();
    auto vertex = tpc.fRPs.front();
    auto directionToSort = line.GetPoint() - vertex;
    lightCl.SortAlongDir(directionToSort);
    lightCl.ReFit();
    auto pos = lightCl.GetVoxels().back().GetPosition();
    return {pos.X(), pos.Y()};
}

// Overlay: "low" (rojo, línea) vs "rest" (negro, puntos), ambos a área unidad
void DrawOverlay(TH1* hLow, TH1* hRest)
{
    hLow->SetStats(false);
    hRest->SetStats(false);
    hLow->SetLineColor(kRed);
    hLow->SetLineWidth(2);
    hRest->SetLineColor(kBlack);
    hRest->SetMarkerColor(kBlack);
    hRest->SetMarkerStyle(20);
    hRest->SetMarkerSize(0.6);

    if(hLow->Integral() > 0)
        hLow->Scale(1. / hLow->Integral());
    if(hRest->Integral() > 0)
        hRest->Scale(1. / hRest->Integral());

    hLow->SetMaximum(1.2 * std::max(hLow->GetMaximum(), hRest->GetMaximum()));
    hLow->GetYaxis()->SetTitle("Normalised counts");

    hLow->DrawClone("HIST");
    hRest->DrawClone("E SAME");

    auto leg = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg->AddEntry(hLow, Form("RawTL < %.0f mm", rawTLCut), "l");
    leg->AddEntry(hRest, Form("RawTL #geq %.0f mm", rawTLCut), "lep");
    leg->Draw();
}

void EventsLowRawTLElastic()
{
    PrettyStyle(true, false);
    // Reacción: 7Li(d,p)8Li, solo L1
    std::string beam {"7Li"};
    std::string targetExp {"d"};
    std::string lightExp {"p"};

    TString fileNameExp {TString::Format("../../PostAnalysis/Outputs/tree_ex_F_%s_%s_%s.root", beam.c_str(),
                                         targetExp.c_str(), lightExp.c_str())};

    ROOT::EnableImplicitMT();

    ROOT::RDataFrame dfExpRaw("Final_Tree", fileNameExp);

    // Columnas originales del árbol (para guardarlas en el Snapshot)
    auto outColumns = dfExpRaw.GetColumnNames();

    // ---------------- Selección L1 + variables ----------------
    auto dfL1 = dfExpRaw.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"})
                    .Filter([](double ex) { return ex < 0.5; }, {"Ex"});

    auto dfExp =
        dfL1.Define("thetaLabDeg", [](ActRoot::MergerData& m) { return static_cast<double>(m.fThetaLight); },
                    {"MergerData"})
            .Define("phiLabDeg", [](ActRoot::MergerData& m) { return static_cast<double>(m.fPhiLight); },
                    {"MergerData"})
            .Define("RPx", [](ActRoot::MergerData& m) { return static_cast<double>(m.fRP.X()); }, {"MergerData"})
            .Define("RPy", [](ActRoot::MergerData& m) { return static_cast<double>(m.fRP.Y()); }, {"MergerData"})
            .Define("RPz", [](ActRoot::MergerData& m) { return static_cast<double>(m.fRP.Z()); }, {"MergerData"})
            .Define("TLExp", [](ActRoot::MergerData& m) { return static_cast<double>(m.fLight.fTL); }, {"MergerData"})
            .Define("RawTLExp", [](ActRoot::MergerData& m) { return static_cast<double>(m.fLight.fRawTL); },
                    {"MergerData"})
            .Define("QtotalLight", [](ActRoot::MergerData& m) { return static_cast<double>(m.fLight.fQtotal); },
                    {"MergerData"})
            .Define("nVoxLight", [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                    { return static_cast<int>(tpc.fClusters[m.fLightIdx].GetVoxels().size()); },
                    {"MergerData", "TPCData"})
            .Define("LastPosX", [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                    { return padSize * GetLastVoxelXY(m, tpc).first; }, {"MergerData", "TPCData"})
            .Define("LastPosY", [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                    { return padSize * GetLastVoxelXY(m, tpc).second; }, {"MergerData", "TPCData"})
            .Define("TLoverRaw", [](double tl, double raw) { return raw > 0 ? tl / raw : -1.; }, {"TLExp", "RawTLExp"});

    // ---------------- Muestras: low vs rest ----------------
    auto dfLow = dfExp.Filter([](double raw) { return raw < rawTLCut; }, {"RawTLExp"});
    auto dfRest = dfExp.Filter([](double raw) { return raw >= rawTLCut; }, {"RawTLExp"});

    // ---------------- Guardado de los eventos low  ----------------
    std::ofstream out("./Outputs/events_Low_RawTL_elastic.dat");
    dfLow.Foreach(
        [&](ActRoot::MergerData& m)
        {
            // Guardar solo los eventos low
            m.Stream(out);
        },
        {"MergerData"});

    // ---------------- Histogramas ----------------
    // RawTL completo (para ver el corte)
    auto hRawAll = dfExp.Histo1D({"hRawAll", "RawTL (L1);RawTL (mm);Counts", 100, 0, 200}, "RawTLExp");

    // Variables 1D: low vs rest
    auto hThetaLow = dfLow.Histo1D(
        {"hThetaLow", "#theta_{Lab};#theta_{Lab} (deg);Counts", thetaBins, thetaMin, thetaMax}, "thetaLabDeg");
    auto hThetaRest = dfRest.Histo1D({"hThetaRest", "", thetaBins, thetaMin, thetaMax}, "thetaLabDeg");

    auto hPhiLow = dfLow.Histo1D({"hPhiLow", "#phi;#phi (deg);Counts", 72, -180, 180}, "phiLabDeg");
    auto hPhiRest = dfRest.Histo1D({"hPhiRest", "", 72, -180, 180}, "phiLabDeg");

    auto hRPxLow = dfLow.Histo1D({"hRPxLow", "RP.X;RP.X (mm);Counts", 128, 0, detMax}, "RPx");
    auto hRPxRest = dfRest.Histo1D({"hRPxRest", "", 128, 0, detMax}, "RPx");
    auto hRPyLow = dfLow.Histo1D({"hRPyLow", "RP.Y;RP.Y (mm);Counts", 128, 0, detMax}, "RPy");
    auto hRPyRest = dfRest.Histo1D({"hRPyRest", "", 128, 0, detMax}, "RPy");
    auto hRPzLow = dfLow.Histo1D({"hRPzLow", "RP.Z;RP.Z (mm);Counts", 128, 0, detMax}, "RPz");
    auto hRPzRest = dfRest.Histo1D({"hRPzRest", "", 128, 0, detMax}, "RPz");

    auto hTLLow = dfLow.Histo1D({"hTLLow", "TL;TL (mm);Counts", 100, 0, 200}, "TLExp");
    auto hTLRest = dfRest.Histo1D({"hTLRest", "", 100, 0, 200}, "TLExp");

    auto hQLow = dfLow.Histo1D({"hQLow", "Q_{total} light;Q_{total};Counts", 100, 0, qMax}, "QtotalLight");
    auto hQRest = dfRest.Histo1D({"hQRest", "", 100, 0, qMax}, "QtotalLight");

    auto hNVoxLow = dfLow.Histo1D({"hNVoxLow", "N voxels light;N voxels;Counts", 100, 0, 200}, "nVoxLight");
    auto hNVoxRest = dfRest.Histo1D({"hNVoxRest", "", 100, 0, 200}, "nVoxLight");

    auto hEVLow = dfLow.Histo1D({"hEVLow", "E_{Vertex};E_{Vertex} (MeV);Counts", 120, 0, eVertexMax}, "EVertex");
    auto hEVRest = dfRest.Histo1D({"hEVRest", "", 120, 0, eVertexMax}, "EVertex");

    auto hExLow = dfLow.Histo1D({"hExLow", "E_{x};E_{x} (MeV);Counts", 100, -5, 5}, "Ex");
    auto hExRest = dfRest.Histo1D({"hExRest", "", 100, -5, 5}, "Ex");

    // 2D solo de la muestra low (y theta vs RawTL de todo)
    auto hRawThetaAll = dfExp.Histo2D({"hRawThetaAll", "RawTL vs #theta_{Lab} (all L1);#theta_{Lab} (deg);RawTL (mm)",
                                       thetaBins, thetaMin, thetaMax, 100, 0, 200},
                                      "thetaLabDeg", "RawTLExp");
    auto hRawQAll =
        dfExp.Histo2D({"hRawQAll", "RawTL vs Q_{total} (all L1);RawTL (mm);Q_{total}", 100, 0, 200, 100, 0, qMax},
                      "RawTLExp", "QtotalLight");
    auto hRawTLAll2D = dfExp.Histo2D(
        {"hRawTLAll2D", "RawTL vs TL (all L1);TL (mm);RawTL (mm)", 100, 0, 200, 100, 0, 200}, "TLExp", "RawTLExp");
    auto hLastLow = dfLow.Histo2D(
        {"hLastLow", "Last voxel (low RawTL);LastPosX (mm);LastPosY (mm)", 128, 0, detMax, 128, 0, detMax}, "LastPosX",
        "LastPosY");
    auto hRPxyLow =
        dfLow.Histo2D({"hRPxyLow", "RP (low RawTL);RP.X (mm);RP.Y (mm)", 128, 0, detMax, 128, 0, detMax}, "RPx", "RPy");
    auto hKinLow = dfLow.Histo2D({"hKinLow", "Kinematics (low RawTL);#theta_{Lab} (deg);E_{Vertex} (MeV)", thetaBins,
                                  thetaMin, thetaMax, 120, 0, eVertexMax},
                                 "thetaLabDeg", "EVertex");

    // Contadores
    auto nLow = dfLow.Count();
    auto nAll = dfExp.Count();

    // Esto dispara el bucle sobre los datos (todo lo anterior es lazy)
    std::cout << "Eventos L1 totales : " << nAll.GetValue() << '\n';
    std::cout << "Eventos RawTL < " << rawTLCut << " : " << nLow.GetValue() << '\n';

    // ---------------- Plots ----------------
    auto c1 = new TCanvas("c1", "RawTL bajo: cinematica y angulos", 1200, 800);
    c1->Divide(2, 2);
    c1->cd(1);
    hRawAll->DrawClone("HIST");
    auto* line = new TLine(rawTLCut, 0, rawTLCut, hRawAll->GetMaximum());
    line->SetLineColor(kRed);
    line->SetLineStyle(2);
    line->Draw();
    c1->cd(2);
    DrawOverlay(hThetaLow.GetPtr(), hThetaRest.GetPtr());
    c1->cd(3);
    DrawOverlay(hPhiLow.GetPtr(), hPhiRest.GetPtr());
    c1->cd(4);
    DrawOverlay(hExLow.GetPtr(), hExRest.GetPtr());

    auto c2 = new TCanvas("c2", "RawTL bajo: punto de reaccion", 1200, 500);
    c2->Divide(3, 1);
    c2->cd(1);
    DrawOverlay(hRPxLow.GetPtr(), hRPxRest.GetPtr());
    c2->cd(2);
    DrawOverlay(hRPyLow.GetPtr(), hRPyRest.GetPtr());
    c2->cd(3);
    DrawOverlay(hRPzLow.GetPtr(), hRPzRest.GetPtr());

    auto c3 = new TCanvas("c3", "RawTL bajo: traza", 1200, 800);
    c3->Divide(2, 2);
    c3->cd(1);
    DrawOverlay(hTLLow.GetPtr(), hTLRest.GetPtr());
    c3->cd(2);
    DrawOverlay(hQLow.GetPtr(), hQRest.GetPtr());
    c3->cd(3);
    DrawOverlay(hNVoxLow.GetPtr(), hNVoxRest.GetPtr());
    c3->cd(4);
    DrawOverlay(hEVLow.GetPtr(), hEVRest.GetPtr());

    auto c4 = new TCanvas("c4", "RawTL: correlaciones (todos L1)", 1200, 500);
    c4->Divide(3, 1);
    c4->cd(1);
    hRawThetaAll->DrawClone("colz");
    c4->cd(2);
    hRawQAll->DrawClone("colz");
    c4->cd(3);
    hRawTLAll2D->DrawClone("colz");

    auto c5 = new TCanvas("c5", "RawTL bajo: posiciones 2D", 1200, 500);
    c5->Divide(3, 1);
    c5->cd(1);
    hLastLow->DrawClone("colz");
    c5->cd(2);
    hRPxyLow->DrawClone("colz");
    c5->cd(3);
    hKinLow->DrawClone("colz");
}