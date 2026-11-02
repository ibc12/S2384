// =============================================================================
//  Compare_automatic_calculations_R_PartialWaves.C
//
//  Dibuja las distribuciones angulares de DWUCK4 (xs_dw4.n: dos columnas
//  "angulo  seccion_eficaz", sin cabecera) producidas por MakeDwuckConv.C.
//
//  Se le pasan LISTAS de Ex, de Rmax y de Lmax y la macro genera, para cada Ex,
//  un canvas con todas las combinaciones (Rmax x Lmax) y lo guarda en png.
//
//     color        -> Lmax
//     tipo de línea -> Rmax
//
//  Rutas que lee (las mismas que escribe MakeDwuckConv.C):
//     <kDirPrefix><tag>/R<R>_L<L>/xs_dw4.n        tag = "gs"  o  "Ex3.000", ...
//
//  Uso:
//     root -l 'Compare_automatic_calculations_R_PartialWaves.C+'
//     root -l 'Compare_automatic_calculations_R_PartialWaves.C+({0.,3.},{30,35,45},{40,50})'
//     root -l 'Compare_automatic_calculations_R_PartialWaves.C+({0.},{20,25,30},{30,40,50,60},false)'  // Y lineal
//
//  Si falta algún archivo (p.ej. esa combinación no se calculó) se avisa y se omite.
// =============================================================================

#include "Rtypes.h"

#include <TAxis.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace CompDW
{

// ----------------------------- CONFIGURACIÓN ---------------------------------
// Carpeta de cada Ex = kDirPrefix + tag   (como crea MakeDwuckConv.C)
const char* const kDirPrefix = "./conv_0p3-2_2+_";
const char* const kXsFile = "xs_dw4.n";
const char* const kFigDir = "./Figures";

// Texto del eje Y / X
const char* const kXTitle = "#theta_{cm} (deg)";
const char* const kYTitle = "d#sigma/d#Omega (mb/sr)";

// ------------------------------ utilidades -----------------------------------
struct Entry
{
    TString path;
    TString label;
    Color_t color;
    Style_t style;
};

// Mismo criterio que MakeDwuckConv.C: "gs" para Ex=0, "Ex3.000" para Ex=3
inline TString Tag(double ex)
{
    if(std::fabs(ex) < 1e-9)
        return "gs";
    return TString(Form("Ex%.3f", ex));
}

inline TString ExTitle(double ex)
{
    if(std::fabs(ex) < 1e-9)
        return "g.s.";
    return TString(Form("E_{x} = %g MeV", ex));
}

inline Color_t ColorL(size_t i)
{
    static const Color_t pal[] = {kBlack,       kRed + 1,    kBlue + 1, kGreen + 2,
                                  kMagenta + 1, kOrange + 7, kCyan + 2, kViolet + 1};
    return pal[i % (sizeof(pal) / sizeof(pal[0]))];
}

inline Style_t StyleR(size_t i)
{
    static const Style_t sty[] = {1, 2, 7, 3, 9, 5};
    return sty[i % (sizeof(sty) / sizeof(sty[0]))];
}

// Lee "angulo xs" (acepta exponentes Fortran 1.0D-03). Líneas no numéricas se ignoran.
// Con logY se descartan los puntos con xs <= 0.
inline bool ReadXs(const TString& path, std::vector<double>& x, std::vector<double>& y, bool logY)
{
    std::ifstream in(path.Data());
    if(!in)
        return false;
    std::string line;
    while(std::getline(in, line))
    {
        for(size_t i = 0; i < line.size(); ++i)
            if(line[i] == 'D' || line[i] == 'd')
                line[i] = 'E';
        std::istringstream ss(line);
        double a, b;
        if(!(ss >> a >> b))
            continue;
        if(logY && b <= 0.)
            continue;
        x.push_back(a);
        y.push_back(b);
    }
    return !x.empty();
}

// ------------------------- canvas con una lista de curvas ---------------------
inline TCanvas* DrawEntries(const std::vector<Entry>& entries, const TString& title, const TString& name,
                            const TString& saveAs, bool logY)
{
    gStyle->SetOptStat(0);

    auto* c = new TCanvas(name, name, 900, 650);
    c->SetLeftMargin(0.13);
    c->SetGrid();
    if(logY)
        c->SetLogy();

    auto* mg = new TMultiGraph(Form("mg_%s", name.Data()), ""); // es el propietario de los TGraph

    std::vector<TGraph*> gs;
    for(size_t i = 0; i < entries.size(); ++i)
    {
        std::vector<double> x, y;
        if(!ReadXs(entries[i].path, x, y, logY))
        {
            std::cerr << "[AVISO] sin datos / no existe: " << entries[i].path << " (se omite)" << std::endl;
            continue;
        }
        auto* g = new TGraph(static_cast<int>(x.size()), x.data(), y.data());
        g->SetName(Form("g_%s_%d", name.Data(), static_cast<int>(i)));
        g->SetTitle(entries[i].label);
        g->SetLineColor(entries[i].color);
        g->SetLineStyle(entries[i].style);
        g->SetLineWidth(2);
        g->SetMarkerStyle(0);
        g->SetFillStyle(0);
        mg->Add(g, "L");
        gs.push_back(g);
    }

    if(gs.empty())
    {
        std::cerr << "[ERROR] " << name << ": ningún archivo pudo leerse, no hay nada que dibujar." << std::endl;
        return c;
    }

    // Leyenda: 2 columnas si hay muchas curvas
    const int n = static_cast<int>(gs.size());
    const int ncol = (n > 8) ? 2 : 1;
    const int rows = (n + ncol - 1) / ncol;
    auto* leg = new TLegend(0.88 - 0.27 * ncol, 0.88 - 0.045 * rows, 0.88, 0.88);
    leg->SetNColumns(ncol);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    for(size_t i = 0; i < gs.size(); ++i)
        leg->AddEntry(gs[i], gs[i]->GetTitle(), "L");

    mg->SetTitle(Form("%s;%s;%s", title.Data(), kXTitle, kYTitle));
    mg->Draw("A");
    mg->GetXaxis()->SetTitleOffset(1.1);
    mg->GetYaxis()->SetTitleOffset(1.3);
    leg->Draw();

    c->Modified();
    c->Update();
    if(saveAs.Length() > 0)
        c->SaveAs(saveAs);
    return c;
}

} // namespace CompDW

// -----------------------------------------------------------------------------
//  Punto de entrada: listas de Ex, Rmax y Lmax
// -----------------------------------------------------------------------------
void Compare_automatic_calculations_R_PartialWaves(const std::vector<double>& exList = {0.0, 3.0},
                                                   const std::vector<int>& rList = {30, 35, 40, 45, 50, 55, 60},
                                                   const std::vector<int>& lList = {40, 45, 50, 55, 60, 65, 70}, bool logY = true)
{
    using namespace CompDW;

    gSystem->mkdir(kFigDir, true);

    for(size_t k = 0; k < exList.size(); ++k)
    {
        const double ex = exList[k];
        const TString tag = Tag(ex);
        const TString dir = TString(kDirPrefix) + tag;

        // L en el bucle externo, R en el interno -> la leyenda agrupa por color (Lmax)
        std::vector<Entry> entries;
        for(size_t j = 0; j < lList.size(); ++j)
            for(size_t i = 0; i < rList.size(); ++i)
            {
                Entry e;
                e.path = Form("%s/R%d_L%d/%s", dir.Data(), rList[i], lList[j], kXsFile);
                e.label = Form("R_{max}=%d fm, L_{max}=%d", rList[i], lList[j]);
                e.color = ColorL(j);
                e.style = StyleR(i);
                entries.push_back(e);
            }

        DrawEntries(entries, ExTitle(ex), TString("cConv_") + tag, TString(kFigDir) + "/compare_" + tag + ".png", logY);
    }
}