// =============================================================================
//  Compare_calc.C
//
//  Macro de ROOT para comparar gráficamente distribuciones angulares
//  (ángulo - sección eficaz) calculadas con TWOFNR y FRESCO (y DWUCK4).
//
//  Formatos de entrada:
//    DWUCK4 : dos columnas "angulo  seccion_eficaz", sin cabecera.
//    FRESCO : mismo formato, con 10 líneas de título al principio
//             y una línea "END" al final.
//    TWOFNR : tres columnas "angulo  seccion_eficaz  (tercera col.)", sin
//             cabecera. Solo se leen las dos primeras columnas.
//
//  Uso:
//    root -l Compare_calc.C
//    root -l 'Compare_calc.C+'     // compilada
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

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace CompCalc
{

// -----------------------------------------------------------------------------
//  Tipos
// -----------------------------------------------------------------------------

enum class Format
{
    DWUCK4,
    FRESCO,
    TWOFNR
};

inline const char* CodeName(Format f)
{
    switch(f)
    {
    case Format::FRESCO: return "FRESCO";
    case Format::TWOFNR: return "TWOFNR";
    default: return "DWUCK4";
    }
}

// Descripción de un cálculo (un archivo) a dibujar
struct Calc
{
    TString path;      // ruta al archivo
    TString label;     // texto de la leyenda
    Format format;     // DWUCK4, FRESCO o TWOFNR
    Color_t color;     // -1 => se asigna automáticamente
    Style_t lineStyle; // 1 = continua, 2 = discontinua, ...
    Width_t lineWidth;
    TString tag; // texto que distingue este cálculo de otros del mismo código

    int skipLines = -1; // -1 => valor por defecto del formato
    int xCol = 0;       // columna del ángulo
    int yCol = 1;       // columna de la sección eficaz

    Calc(const TString& p, const TString& l, Format f, Color_t c = -1, Style_t ls = 1, Width_t lw = 2,
         const TString& t = "")
        : path(p),
          label(l),
          format(f),
          color(c),
          lineStyle(ls),
          lineWidth(lw),
          tag(t)
    {
    }

    // Setters encadenables
    Calc& Tag(const TString& t)
    {
        tag = t;
        return *this;
    }
    Calc& Skip(int n)
    {
        skipLines = n;
        return *this;
    }
    Calc& Cols(int ix, int iy)
    {
        xCol = ix;
        yCol = iy;
        return *this;
    }

    // Texto final de la leyenda:
    //   label vacío  -> nombre del código
    //   tag no vacío -> se añade entre paréntesis: "FRESCO (pot. Koning)"
    TString LegendText() const
    {
        TString s = label;
        if(s.IsNull())
            s = CodeName(format);
        if(!tag.IsNull())
            s += " (" + tag + ")";
        return s;
    }
};

// Opciones globales de la figura
struct PlotOptions
{
    TString title = "";
    TString xTitle = "#theta_{cm} (deg)";
    TString yTitle = "d#sigma/d#Omega (mb/sr)";
    bool logY = true;
    TString canvasName = "cComp";
    int width = 900;
    int height = 650;
    TString saveAs = "";
};

// Parámetros de cada formato
constexpr int kFrescoHeaderLines = 10;

inline int DefaultHeaderLines(Format f)
{
    return f == Format::FRESCO ? kFrescoHeaderLines : 0; // DWUCK4 y TWOFNR: sin cabecera
}
inline bool HasEndMarker(Format f)
{
    return f == Format::FRESCO; // DWUCK4 y TWOFNR: sin "END"
}

// -----------------------------------------------------------------------------
//  Utilidades de texto
// -----------------------------------------------------------------------------

inline std::string Trim(const std::string& s)
{
    const char* ws = " \t\r\n";
    const auto a = s.find_first_not_of(ws);
    if(a == std::string::npos)
        return "";
    const auto b = s.find_last_not_of(ws);
    return s.substr(a, b - a + 1);
}

inline bool IsEndMarker(const std::string& line)
{
    std::string t = Trim(line);
    std::transform(t.begin(), t.end(), t.begin(), ::toupper);
    return t == "END";
}

// Lee todos los números de una línea (acepta exponentes Fortran 1.0D-03).
// Devuelve true si la línea tiene al menos 'nNeeded' números válidos.
inline bool ParseNumbers(std::string line, std::vector<double>& vals, size_t nNeeded)
{
    for(auto& c : line)
        if(c == 'D' || c == 'd')
            c = 'E';
    std::istringstream ss(line);
    double v;
    vals.clear();
    while(ss >> v)
        vals.push_back(v);
    return vals.size() >= nNeeded;
}

// -----------------------------------------------------------------------------
//  Lectura de archivos
// -----------------------------------------------------------------------------

// Devuelve true si se ha leído al menos un punto.
inline bool ReadData(const Calc& calc, std::vector<double>& x, std::vector<double>& y)
{
    std::ifstream in(calc.path.Data());
    if(!in)
    {
        std::cerr << "[ERROR] No se puede abrir: " << calc.path << std::endl;
        return false;
    }

    const int skip = (calc.skipLines >= 0) ? calc.skipLines : DefaultHeaderLines(calc.format);
    const bool useEnd = HasEndMarker(calc.format);
    const size_t nNeeded = static_cast<size_t>(std::max(calc.xCol, calc.yCol)) + 1;

    std::string line;
    std::vector<double> vals;
    int nLine = 0, nBad = 0;

    while(std::getline(in, line))
    {
        ++nLine;
        if(nLine <= skip)
            continue; // cabecera
        if(useEnd && IsEndMarker(line))
            break; // fin de datos
        if(Trim(line).empty())
            continue; // línea vacía

        if(!ParseNumbers(line, vals, nNeeded))
        {
            ++nBad;
            continue;
        }
        x.push_back(vals[calc.xCol]);
        y.push_back(vals[calc.yCol]);
    }

    if(nBad > 0)
        std::cerr << "[AVISO] " << calc.path << ": " << nBad << " línea(s) no numérica(s) ignorada(s)." << std::endl;
    if(x.empty())
    {
        std::cerr << "[ERROR] Sin datos válidos en: " << calc.path << std::endl;
        return false;
    }
    return true;
}

// -----------------------------------------------------------------------------
//  Construcción del TGraph
// -----------------------------------------------------------------------------

// Si logY == true, se descartan los puntos con y <= 0 (no representables).
inline TGraph* MakeGraph(const Calc& calc, Color_t color, bool logY)
{
    std::vector<double> x, y;
    if(!ReadData(calc, x, y))
        return nullptr;

    if(logY)
    {
        std::vector<double> xf, yf;
        for(size_t i = 0; i < x.size(); ++i)
            if(y[i] > 0)
            {
                xf.push_back(x[i]);
                yf.push_back(y[i]);
            }
        if(xf.size() != x.size())
            std::cerr << "[AVISO] " << calc.path << ": " << (x.size() - xf.size())
                      << " punto(s) con y<=0 descartado(s) por escala log." << std::endl;
        x.swap(xf);
        y.swap(yf);
        if(x.empty())
            return nullptr;
    }

    auto* g = new TGraph(static_cast<int>(x.size()), x.data(), y.data());
    g->SetName(Form("g_%s", gSystem->BaseName(calc.path)));
    g->SetTitle(calc.LegendText());
    g->SetLineColor(color);
    g->SetLineStyle(calc.lineStyle);
    g->SetLineWidth(calc.lineWidth);
    g->SetMarkerStyle(0);
    g->SetFillStyle(0);
    return g;
}

// -----------------------------------------------------------------------------
//  Figura completa
// -----------------------------------------------------------------------------

inline Color_t AutoColor(size_t i)
{
    static const Color_t palette[] = {kBlack,      kRed + 1,  kBlue + 1,   kGreen + 2, kMagenta + 1,
                                      kOrange + 7, kCyan + 2, kViolet + 1, kGray + 2,  kYellow + 2};
    return palette[i % (sizeof(palette) / sizeof(palette[0]))];
}

// Dibuja todos los cálculos en un único canvas. Devuelve el canvas.
inline TCanvas* DrawComparison(const std::vector<Calc>& calcs, const PlotOptions& opt = PlotOptions())
{
    gStyle->SetOptStat(0);

    auto* c = new TCanvas(opt.canvasName, opt.canvasName, opt.width, opt.height);
    c->SetLeftMargin(0.13);
    c->SetGrid();
    if(opt.logY)
        c->SetLogy();

    // El TMultiGraph es el propietario de los TGraph que se le añaden
    auto* mg = new TMultiGraph(Form("mg_%s", opt.canvasName.Data()), "");

    // Ancho de la leyenda según el texto más largo (entre 25% y 70% del canvas)
    size_t maxLen = 0;
    for(const auto& cc : calcs)
        maxLen = std::max<size_t>(maxLen, cc.LegendText().Length());
    const double legW = std::min(0.70, std::max(0.25, 0.10 + 0.012 * maxLen));
    auto* leg = new TLegend(0.88 - legW, 0.88 - 0.05 * calcs.size(), 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    int nOk = 0;
    for(size_t i = 0; i < calcs.size(); ++i)
    {
        const Color_t col = (calcs[i].color >= 0) ? calcs[i].color : AutoColor(i);
        TGraph* g = MakeGraph(calcs[i], col, opt.logY);
        if(!g)
            continue;
        mg->Add(g, "L");
        leg->AddEntry(g, calcs[i].LegendText(), "L");
        ++nOk;
    }

    if(nOk == 0)
    {
        std::cerr << "[ERROR] Ningún archivo pudo leerse; no hay nada que dibujar." << std::endl;
        return c;
    }

    mg->SetTitle(Form("%s;%s;%s", opt.title.Data(), opt.xTitle.Data(), opt.yTitle.Data()));
    mg->Draw("A");
    mg->GetXaxis()->SetTitleOffset(1.1);
    mg->GetYaxis()->SetTitleOffset(1.3);
    leg->Draw();

    c->Modified();
    c->Update();
    if(opt.saveAs.Length() > 0)
        c->SaveAs(opt.saveAs);
    return c;
}

} // namespace CompCalc

// -----------------------------------------------------------------------------
//  Punto de entrada: EDITA AQUÍ tus archivos (TWOFNR vs FRESCO)
//  TWOFNR: línea continua.  FRESCO: línea discontinua.
// -----------------------------------------------------------------------------
void Compare_calc_FRESCO_TWOFNR()
{
    using namespace CompCalc;

    // --- Estado fundamental ---
    std::vector<Calc> gs = {
        Calc("./gs_ADWA_KD_twofnr_JohnSoper_LEA(zerorange)/21.gs", "TWOFNR", Format::TWOFNR, kBlack).Tag("ADKD-ZR-LEA"),
        Calc("./gs_ADWA_KD_twofnr_JohnTandy_LEA(finiterange)/21.gs", "TWOFNR", Format::TWOFNR, kPink)
            .Tag("ADKD-FR-LEA"),
        Calc("./gs_ADWA_KD_twofnr_JohnSoper_LEA(zerorange)/fort.202", "FRESCO", Format::FRESCO, kBlack, 2)
            .Tag("ADKD-ZR-LEA (in twofnr, ZR in fresco)"),
        Calc("./gs_ADWA_KD_twofnr_JohnTandy_LEA(finiterange)/fort.202", "FRESCO", Format::FRESCO, kPink, 2)
            .Tag("ADKD-FR-LEA (in twofnr, ZR in fresco)"),
    };
    PlotOptions optGS;
    optGS.title = "g.s.";
    optGS.canvasName = "cGS";
    optGS.saveAs = "./Figures/compare_twofnr_gs.png";
    DrawComparison(gs, optGS);

    // --- Estado fundamental SO contribution---
    std::vector<Calc> gs_SO = {
        // --- Zero range (línea discontinua) ---
        Calc("./gs_ADWA_KD_twofnr_JohnSoper_LEA(zerorange)/fort.202", "FRESCO", Format::FRESCO, kBlack, 2)
            .Tag("ADKD-ZR-LEA, (in twofnr, ZR in fresco), SO completo "),
        Calc("./gs_ADWA_KD_twofnr_JohnSoper_LEA(zerorange)_noSO_entrance/fort.202", "FRESCO", Format::FRESCO, kRed + 1,
             2)
            .Tag("ADKD-ZR-LEA, (in twofnr, ZR in fresco),noSO entrance "),
        Calc("./gs_ADWA_KD_twofnr_JohnSoper_LEA(zerorange)_noSO_exit/fort.202", "FRESCO", Format::FRESCO, kBlue + 1, 2)
            .Tag("ADKD-ZR-LEA,(in twofnr, ZR in fresco), noSO exit"),
        Calc("./gs_ADWA_KD_twofnr_JohnSoper_LEA(zerorange)_noSO_both/fort.202", "FRESCO", Format::FRESCO, kGreen + 2, 2)
            .Tag("ADKD-ZR-LEA,(in twofnr, ZR in fresco), noSO both"),
        // --- Finite range (línea continua) ---
        Calc("./gs_ADWA_KD_twofnr_JohnTandy_LEA(finiterange)/fort.202", "FRESCO", Format::FRESCO, kBlack, 1)
            .Tag("ADKD-FR-LEA, (in twofnr, ZR in fresco), SO completo "),
        Calc("./gs_ADWA_KD_twofnr_JohnTandy_LEA(finiterange)_noSO_entrance/fort.202", "FRESCO", Format::FRESCO,
             kRed + 1, 1)
            .Tag("ADKD-FR-LEA, (in twofnr, ZR in fresco), noSO entrance "),
        Calc("./gs_ADWA_KD_twofnr_JohnTandy_LEA(finiterange)_noSO_exit/fort.202", "FRESCO", Format::FRESCO, kBlue + 1,
             1)
            .Tag("ADKD-FR-LEA, (in twofnr, ZR in fresco), noSO exit"),
        Calc("./gs_ADWA_KD_twofnr_JohnTandy_LEA(finiterange)_noSO_both/fort.202", "FRESCO", Format::FRESCO, kGreen + 2,
             1)
            .Tag("ADKD-FR-LEA,(in twofnr, ZR in fresco), noSO both"),
    };
    PlotOptions optGS_SO;
    optGS_SO.title = "g.s. SO contribution";
    optGS_SO.canvasName = "cGS_SO";
    optGS_SO.saveAs = "./Figures/compare_twofnr_gs_SO.png";
    DrawComparison(gs_SO, optGS_SO);

    // --- Primer estado, Ex = 0.981 MeV ---
    // std::vector<Calc> g1 = {
    //     Calc("./TWOFNR_calculations/ADKDWat_noSO/g1/xs_twofnr.dat", "TWOFNR", Format::TWOFNR, kBlack)
    //         .Tag("ADKDWat-no SO"),
    //     Calc("./TWOFNR_calculations/ADKDKD_noSO/g1/xs_twofnr.dat", "TWOFNR", Format::TWOFNR, kPink).Tag("ADKDKD-no
    //     SO"), Calc("./TWOFNR_calculations/ADKDKD/g1/xs_twofnr.dat", "TWOFNR", Format::TWOFNR, kGreen).Tag("ADKDKD"),
    //     Calc("./gs_ADWA_Watson_FRESCO/fort.203", "FRESCO", Format::FRESCO, kBlack, 2).Tag("ADKDWat-no SO"),
    //     Calc("./gs_ADWA_KD_FRESCO_noSO/fort.203", "FRESCO", Format::FRESCO, kPink + 1, 2).Tag("ADKDKD-no SO"),
    //     Calc("./gs_ADWA_KD_FRESCO/fort.203", "FRESCO", Format::FRESCO, kGreen + 1, 2).Tag("ADKDKD"),
    // };
    // PlotOptions optG1;
    // optG1.title = "E_{x} = 0.981 MeV";
    // optG1.canvasName = "cG1";
    // optG1.saveAs = "./Figures/compare_twofnr_0-981MeV.png";
    // DrawComparison(g1, optG1);
}