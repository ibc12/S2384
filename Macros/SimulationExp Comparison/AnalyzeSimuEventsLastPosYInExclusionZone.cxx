#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>

// ---------------- Configuración ----------------
constexpr int yMinExclusionZone = 56;
constexpr int yMaxExclusionZone = 71;
constexpr double padSize = 2.;
constexpr double yMinZoneMm = yMinExclusionZone * padSize;       // 112 mm
constexpr double yMaxZoneMm = (yMaxExclusionZone + 1) * padSize; // 144 mm
constexpr double detMax = 256.;                                  // mm (X e Y: 0-256)

// true -> el Y predicho se invierte (Y_pred -> 256 - Y_pred) para compararlo con LastPosY
constexpr bool flipYPred = true;

inline bool InZone(double y)
{
    return y >= yMinZoneMm && y <= yMaxZoneMm;
}

// Punto de salida de la pesada: {X, Y, motivo}; th4, phi4 en RADIANES
std::array<double, 3> HeavyExit(double vx, double vy, double th4, double ph4)
{
    const double dx = std::cos(th4);
    const double dy = std::sin(th4) * std::sin(ph4);
    if(dx <= 1e-9)
        return {vx, vy, -1.};
    const double inf = std::numeric_limits<double>::infinity();
    const double tX = (detMax - vx) / dx;
    double tY = inf;
    if(dy > 1e-12)
        tY = (detMax - vy) / dy;
    else if(dy < -1e-12)
        tY = (0. - vy) / dy;
    const double t = std::min(tX, tY);
    return {vx + t * dx, vy + t * dy, tY < tX ? 1. : 0.};
}

// Overlay: "all" (negro, puntos) vs "in" (rojo, línea), ambos a área unidad
void DrawOverlay(TH1* hAll, TH1* hIn)
{
    hAll->SetStats(false);
    hIn->SetStats(false);
    hIn->SetLineColor(kRed);
    hIn->SetLineWidth(2);
    hAll->SetLineColor(kBlack);
    hAll->SetMarkerColor(kBlack);
    hAll->SetMarkerStyle(20);
    hAll->SetMarkerSize(0.6);
    if(hAll->Integral() > 0)
        hAll->Scale(1. / hAll->Integral());
    if(hIn->Integral() > 0)
        hIn->Scale(1. / hIn->Integral());
    hAll->SetMaximum(1.2 * std::max(hAll->GetMaximum(), hIn->GetMaximum()));
    hAll->GetYaxis()->SetTitle("Normalised counts");
    hAll->DrawClone("E");
    hIn->DrawClone("HIST SAME");
    auto leg = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg->AddEntry(hAll, "All events", "lep");
    leg->AddEntry(hIn, "Light ends IN zone", "l");
    leg->Draw();
}

// Líneas de la zona de exclusión: horizontales (Y en el eje vertical)
void DrawZoneLinesH(double xmax)
{
    for(double y : {yMinZoneMm, yMaxZoneMm})
    {
        auto l = new TLine(0, y, xmax, y);
        l->SetLineColor(kBlue);
        l->SetLineStyle(2);
        l->Draw();
    }
}

void AnalyzeSimuEventsLastPosYInExclusionZone()
{
    std::string beam {"7Li"};
    std::string target {"2H"};
    std::string light {"2H"};
    double Ex {0.};

    TString fileNameSimu {TString::Format("../../Simulation/Outputs/%s/test_rpCut/"
                                          "%s_%s_TRIUMF_Eex_%.3f_nPS_%d_pPS_%d_L1_240rpCut_paralelization.root",
                                          beam.c_str(), target.c_str(), light.c_str(), Ex, 0, 0)};

    ROOT::EnableImplicitMT();
    gStyle->SetOptStat(0);

    ROOT::RDataFrame dfRaw("SimulationTTree", fileNameSimu);

    // Autodetección de unidades (grados si el máximo > 2pi)
    auto toRadFor = [&](const char* col)
    {
        const double mx = *dfRaw.Max(col);
        const bool deg = mx > 2. * M_PI + 1e-6;
        std::cout << col << " max = " << mx << " -> " << (deg ? "grados" : "radianes") << '\n';
        return deg ? M_PI / 180. : 1.;
    };
    const double toRadTh4 = toRadFor("theta4Lab");
    const double toRadPh4 = toRadFor("phi4Lab");
    const double toRadTh3 = toRadFor("theta3Lab");
    const double toRadPh3 = toRadFor("phi3CM");

    // ---------------- Variables ----------------
    auto df = dfRaw
                  .Define("heavyExit", [=](double vx, double vy, double th4, double ph4)
                          { return HeavyExit(vx, vy, th4 * toRadTh4, ph4 * toRadPh4); },
                          {"vertexX", "vertexY", "theta4Lab", "phi4Lab"})
                  .Define("heavyExitX", [](const std::array<double, 3>& e) { return e[0]; }, {"heavyExit"})
                  .Define("heavyExitY", [](const std::array<double, 3>& e) { return e[1]; }, {"heavyExit"})
                  .Define("lightInZone", [](double y) { return InZone(y); }, {"LastPosY"})
                  .Define("theta4LabDeg", [=](double t) { return t * toRadTh4 * 180. / M_PI; }, {"theta4Lab"})
                  // Punto predicho: vértice + TL * dirección de la luz
                  .Define("predX", [=](double vx, double L, double t) { return vx + L * std::cos(t * toRadTh3); },
                          {"vertexX", "TL", "theta3Lab"})
                  .Define("predY",
                          [=](double vy, double L, double t, double p)
                          {
                              const double y = vy + L * std::sin(t * toRadTh3) * std::sin(p * toRadPh3);
                              return flipYPred ? detMax - y : y;
                          },
                          {"vertexY", "TL", "theta3Lab", "phi3CM"});

    auto dfIn = df.Filter("lightInZone");

    // ---------------- Histogramas ----------------
    // 1) Último punto X vs Y: luz y pesada, todos / IN
    auto hLightAll =
        df.Histo2D({"hLightAll", "Light last point (all);LastPosX (mm);LastPosY (mm)", 128, 0, detMax, 128, 0, detMax},
                   "LastPosX", "LastPosY");
    auto hLightIn = dfIn.Histo2D(
        {"hLightIn", "Light last point (IN zone);LastPosX (mm);LastPosY (mm)", 128, 0, detMax, 128, 0, detMax},
        "LastPosX", "LastPosY");
    auto hHeavyAll =
        df.Histo2D({"hHeavyAll", "Heavy last point (all);Heavy X (mm);Heavy Y (mm)", 128, 0, detMax, 128, 0, detMax},
                   "heavyExitX", "heavyExitY");
    auto hHeavyIn = dfIn.Histo2D(
        {"hHeavyIn", "Heavy last point (light IN zone);Heavy X (mm);Heavy Y (mm)", 128, 0, detMax, 128, 0, detMax},
        "heavyExitX", "heavyExitY");

    // 2) LastPosY de la luz vs theta4 (pesada), todos / IN
    auto hYvsTh4All =
        df.Histo2D({"hYvsTh4All", "All;#theta_{4,Lab} (deg);Light LastPosY (mm)", 100, 0, 90, 128, 0, detMax},
                   "theta4LabDeg", "LastPosY");
    auto hYvsTh4In = dfIn.Histo2D(
        {"hYvsTh4In", "Light IN zone;#theta_{4,Lab} (deg);Light LastPosY (mm)", 100, 0, 90, 128, 0, detMax},
        "theta4LabDeg", "LastPosY");

    // 3) Predicho vs real
    auto hPredX =
        df.Histo2D({"hPredX", "LastPosX vs predicted;predX (mm);LastPosX (mm)", 128, 0, detMax, 128, 0, detMax},
                   "predX", "LastPosX");
    auto hPredY =
        df.Histo2D({"hPredY", Form("LastPosY vs predicted%s;predY (mm);LastPosY (mm)", flipYPred ? " (Y flipped)" : ""),
                    128, 0, detMax, 128, 0, detMax},
                   "predY", "LastPosY");

    // 4) nPads: luz, pesada y total; todos / IN
    auto hNLightAll = df.Histo1D({"hNLightAll", "Pads light;N pads;Counts", 150, 0, 150}, "nPadsLight");
    auto hNLightIn = dfIn.Histo1D({"hNLightIn", "", 150, 0, 150}, "nPadsLight");
    auto hNHeavyAll = df.Histo1D({"hNHeavyAll", "Pads heavy;N pads;Counts", 150, 0, 150}, "nPadsHeavy");
    auto hNHeavyIn = dfIn.Histo1D({"hNHeavyIn", "", 150, 0, 150}, "nPadsHeavy");
    auto hNTotAll = df.Histo1D({"hNTotAll", "Pads total;N pads;Counts", 150, 0, 150}, "nPads");
    auto hNTotIn = dfIn.Histo1D({"hNTotIn", "", 150, 0, 150}, "nPads");

    // ---------------- Plots ----------------
    auto c1 = new TCanvas("c1", "Ultimo punto: luz y pesada", 1200, 800);
    c1->Divide(2, 2);
    c1->cd(1);
    hLightAll->DrawClone("colz");
    DrawZoneLinesH(detMax);
    c1->cd(2);
    hLightIn->DrawClone("colz");
    DrawZoneLinesH(detMax);
    c1->cd(3);
    hHeavyAll->DrawClone("colz");
    DrawZoneLinesH(detMax);
    c1->cd(4);
    hHeavyIn->DrawClone("colz");
    DrawZoneLinesH(detMax);

    auto c2 = new TCanvas("c2", "LastPosY luz vs theta pesada", 1200, 500);
    c2->Divide(2, 1);
    c2->cd(1);
    hYvsTh4All->DrawClone("colz");
    DrawZoneLinesH(90);
    c2->cd(2);
    hYvsTh4In->DrawClone("colz");
    DrawZoneLinesH(90);

    auto c3 = new TCanvas("c3", "Predicho vs real", 1200, 500);
    c3->Divide(2, 1);
    c3->cd(1);
    hPredX->DrawClone("colz");
    c3->cd(2);
    hPredY->DrawClone("colz");

    auto c4 = new TCanvas("c4", "nPads luz, pesada, total", 1500, 500);
    c4->Divide(3, 1);
    c4->cd(1);
    DrawOverlay(hNLightAll.GetPtr(), hNLightIn.GetPtr());
    c4->cd(2);
    DrawOverlay(hNHeavyAll.GetPtr(), hNHeavyIn.GetPtr());
    c4->cd(3);
    DrawOverlay(hNTotAll.GetPtr(), hNTotIn.GetPtr());
}