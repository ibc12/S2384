// =============================================================================
//  compare_calc.C
//
//  Macro de ROOT para comparar gráficamente distribuciones angulares
//  (ángulo - sección eficaz) calculadas con DWUCK4 y/o FRESCO.
//
//  Formatos de entrada:
//    DWUCK4 : dos columnas "angulo  seccion_eficaz", sin cabecera.
//    FRESCO : mismo formato, pero con 10 líneas de título al principio
//             y una línea "END" al final.
//
//  Uso:
//    root -l compare_calc.C            // ejecuta compare_calc() con el ejemplo
//    root -l 'compare_calc.C+'         // compilada (más rápido y más estricto)
//
//  Para tus propios archivos, edita la función compare_calc() al final
//  (o llama a CompCalc::DrawComparison(...) desde otra macro).
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
    FRESCO
};

// Descripción de un cálculo (un archivo) a dibujar
struct Calc
{
    TString path;      // ruta al archivo
    TString label;     // texto de la leyenda
    Format format;     // DWUCK4 o FRESCO
    Color_t color;     // -1 => se asigna automáticamente
    Style_t lineStyle; // 1 = continua, 2 = discontinua, ...
    Width_t lineWidth;
    TString tag; // texto que distingue este cálculo de otros
                 // del mismo código (p.ej. "pot. Koning", "r0=1.25")

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

    // Setter encadenable:  Calc(...).Tag("CCBA, sin spin-orbita")
    Calc& Tag(const TString& t)
    {
        tag = t;
        return *this;
    }

    // Texto final de la leyenda:
    //   label vacío  -> nombre del código ("DWUCK4" / "FRESCO")
    //   tag no vacío -> se añade entre paréntesis: "FRESCO (pot. Koning)"
    TString LegendText() const
    {
        TString s = label;
        if(s.IsNull())
            s = (format == Format::FRESCO) ? "FRESCO" : "DWUCK4";
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
    bool logY = true; // escala logarítmica en Y
    TString canvasName = "cComp";
    int width = 900;
    int height = 650;
    TString saveAs = ""; // p.ej. "comparacion.png" o ".pdf"
};

// Parámetros de cada formato
constexpr int kFrescoHeaderLines = 10;

inline int HeaderLines(Format f)
{
    return f == Format::FRESCO ? kFrescoHeaderLines : 0;
}
inline bool HasEndMarker(Format f)
{
    return f == Format::FRESCO;
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

// Lee dos números de una línea. Acepta exponentes estilo Fortran (1.0D-03).
inline bool ParseLine(std::string line, double& x, double& y)
{
    for(auto& c : line)
        if(c == 'D' || c == 'd')
            c = 'E';
    std::istringstream ss(line);
    return static_cast<bool>(ss >> x >> y);
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

    const int skip = HeaderLines(calc.format);
    const bool useEnd = HasEndMarker(calc.format);
    std::string line;
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

        double a, b;
        if(!ParseLine(line, a, b))
        {
            ++nBad;
            continue;
        }
        x.push_back(a);
        y.push_back(b);
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
//  Punto de entrada: EDITA AQUÍ tus archivos
// -----------------------------------------------------------------------------
void Compare_calc_Mario()
{
    using namespace CompCalc;

    // --- All Mario calculations---
    std::vector<Calc> calc_Mario = {
        Calc("../../FranckCalculations/DWUCK4_personal_calculations/good_calculations/v1/xs_dw4.n", "DWUCK4",
             Format::DWUCK4, AutoColor(0), 1)
            .Tag("both SO - my version 3 MeV"),
        Calc("../../FranckCalculations/DWUCK4_examples_Franck/3MeV_Ex_pwave_noOscilation/xs_dw4.n", "DWUCK4",
             Format::DWUCK4, AutoColor(1), 1)
            .Tag("no SO - Franck's version 3 MeV"),
        Calc("./7Lidp8Li/fort.202", "FRESCO", Format::FRESCO, AutoColor(2), 1)
            .Tag("ADWA post - Mario's version near 3 MeV"),
        Calc("./7Lidp8Li_prior/fort.202", "FRESCO", Format::FRESCO, AutoColor(3), 1)
            .Tag("ADWA prior - Mario's version near 3 MeV"),
        Calc("./7Lidp8Li_DWBA/fort.202", "FRESCO", Format::FRESCO, AutoColor(4), 1)
            .Tag("DWBA post - Mario's version near 3 MeV"),
        Calc("./7Lidp8Li_DWBA_prior/fort.202", "FRESCO", Format::FRESCO, AutoColor(5), 1)
            .Tag("DWBA prior - Mario's version near 3 MeV"),
    };
    PlotOptions optMario;
    optMario.title = "All Mario calculations";
    optMario.canvasName = "cMario";
    optMario.saveAs = "./Figures/compare_mario.png";
    DrawComparison(calc_Mario, optMario);


    // --- DWBA prior comparison---
    std::vector<Calc> calc_DWBA_prior = {
        Calc("./7Lidp8Li_DWBA_prior/fort.202", "FRESCO", Format::FRESCO, AutoColor(0), 1)
            .Tag("DWBA prior - Mario's version near 3 MeV"),
        Calc("./personal_calculations/DWBA/prior/DA1p_change/fort.202", "FRESCO", Format::FRESCO, AutoColor(1), 1)
            .Tag("DWBA prior - DA1p change"),
        Calc("./personal_calculations/DWBA/prior/Correct_exitOMP/fort.202", "FRESCO", Format::FRESCO, AutoColor(2), 1)
            .Tag("DWBA prior - Correct exit OMP"),
        Calc("./personal_calculations/DWBA/prior/Correct_remanentTerm/fort.202", "FRESCO", Format::FRESCO, AutoColor(3),
             1)
            .Tag("DWBA prior - Correct remanent term"),
        Calc("./personal_calculations/DWBA/prior/LEA_ZR/fort.202", "FRESCO", Format::FRESCO, AutoColor(4), 1)
            .Tag("DWBA prior - LEA ZR change"),
        Calc("./personal_calculations/DWBA/prior/Change_to_ReidSoftCore/fort.202", "FRESCO", Format::FRESCO,
             AutoColor(5), 1)
            .Tag("DWBA prior - RedSoftCore"),
    };
    PlotOptions optDWBAprior;
    optDWBAprior.title = "DWBA prior comparison";
    optDWBAprior.canvasName = "cDWBAprior";
    optDWBAprior.saveAs = "./Figures/compare_dwba_prior.png";
    DrawComparison(calc_DWBA_prior, optDWBAprior);


    // --- ADWA prior comparison---
    std::vector<Calc> calc_ADWA_prior = {
        Calc("./7Lidp8Li_prior/fort.202", "FRESCO", Format::FRESCO, AutoColor(0), 1)
            .Tag("ADWA prior - Mario's version near 3 MeV"),
        Calc("./personal_calculations/ADWA/prior/myAdiabatic/fort.202", "FRESCO", Format::FRESCO, AutoColor(1), 1)
            .Tag("ADWA prior - my fort.4 from twofnr directly"),
        Calc("./personal_calculations/ADWA/prior/twocols_changeWithMyAdiabatic/fort.202", "FRESCO", Format::FRESCO,
             AutoColor(2), 1)
            .Tag("ADWA prior - 2col change"),
        Calc("./personal_calculations/ADWA/prior/Correct_exitOMP/fort.202", "FRESCO", Format::FRESCO, AutoColor(3), 1)
            .Tag("ADWA prior - Correct exit OMP"),
        Calc("./personal_calculations/ADWA/prior/Correct_remanentTerm/fort.202", "FRESCO", Format::FRESCO, AutoColor(4),
             1)
            .Tag("ADWA prior - Correct remanent term"),
        Calc("./personal_calculations/ADWA/prior/LEA_ZR/fort.202", "FRESCO", Format::FRESCO, AutoColor(5), 1)
            .Tag("ADWA prior - LEA ZR change"),
        Calc("./personal_calculations/ADWA/prior/Change_to_ReidSoftCore/fort.202", "FRESCO", Format::FRESCO,
             AutoColor(6), 1)
            .Tag("ADWA prior - RedSoftCore"),
    };
    PlotOptions optADWAprior;
    optADWAprior.title = "ADWA prior comparison";
    optADWAprior.canvasName = "cADWAprior";
    optADWAprior.saveAs = "./Figures/compare_adwa_prior.png";
    DrawComparison(calc_ADWA_prior, optADWAprior);


    // --- Check elastic calculations ---
    std::vector<Calc> calc_elastic = {
        Calc("../../7Li_dd/Inputs/gsDA1p_corr/fort.201", "FRESCO", Format::FRESCO, AutoColor(0), 1)
            .Tag("elastic - my version - DA1p"),
        Calc("./7Lidp8Li_DWBA_prior/fort.201", "FRESCO", Format::FRESCO, AutoColor(1), 1)
            .Tag("elastic - Mario's DWBA prior"),
        Calc("./7Lidp8Li_prior/fort.201", "FRESCO", Format::FRESCO, AutoColor(2), 1)
            .Tag("elastic - Mario's ADWA prior"),
        Calc("./7Lidp8Li_DWBA_prior/fort.201", "FRESCO", Format::FRESCO, AutoColor(3), 1)
            .Tag("elastic - Mario's DWBA prior with DA1p"),
        Calc("./personal_calculations/DWBA/prior/DA1p_change_sameInputPars(thMax)/fort.201", "FRESCO", Format::FRESCO,
             AutoColor(4), 1)
            .Tag("elastic - Mario's DWBA prior with DA1p and with -180 as thmax (no Rutherfor normalization)"),
        Calc("./personal_calculations/ADWA/prior/Change_sameInputPars(thMax)/fort.201", "FRESCO", Format::FRESCO,
             AutoColor(5), 1)
            .Tag("elastic - Mario's ADWA prior with -180 as thmax (no Rutherfor normalization)"),
    };
    PlotOptions optElastic;
    optElastic.title = "Elastic comparison";
    optElastic.canvasName = "cElastic";
    optElastic.saveAs = "./Figures/compare_elastic.png";
    DrawComparison(calc_elastic, optElastic);

    // --- Total comparison ---
    std::vector<Calc> calc_comparison = {
        Calc("../../FranckCalculations/DWUCK4_personal_calculations/good_calculations/v1/xs_dw4.n", "DWUCK4",
             Format::DWUCK4, AutoColor(0), 1)
            .Tag("both SO - my version 3 MeV"),
        Calc("../../FranckCalculations/DWUCK4_examples_Franck/3MeV_Ex_pwave_noOscilation/xs_dw4.n", "DWUCK4",
             Format::DWUCK4, AutoColor(1), 1)
            .Tag("no SO - Franck's version 3 MeV"),
        Calc("./7Lidp8Li_prior/fort.202", "FRESCO", Format::FRESCO, AutoColor(3), 1)
            .Tag("ADWA prior - Mario's version near 3 MeV"),
        Calc("./7Lidp8Li_DWBA_prior/fort.202", "FRESCO", Format::FRESCO, AutoColor(5), 1)
            .Tag("DWBA prior - Mario's version near 3 MeV"),
        Calc("./personal_calculations/ADWA/prior/ALL_noZR_LEA/fort.202", "FRESCO", Format::FRESCO, AutoColor(6), 1)
            .Tag("ADWA prior - All changes no ZR and LEA"),
        Calc("./personal_calculations/DWBA/prior/ALL_noZR_LEA/fort.202", "FRESCO", Format::FRESCO, AutoColor(7), 1)
            .Tag("DWBA prior - All changes no ZR and LEA"),
        Calc("./personal_calculations/DWBA/prior/ALL/fort.202", "FRESCO", Format::FRESCO, AutoColor(8), 1)
            .Tag("DWBA prior - All changes with ZR and LEA"),
        Calc("./personal_calculations/ADWA/prior/ALL/fort.202", "FRESCO", Format::FRESCO, AutoColor(9), 1)
            .Tag("ADWA prior - All changes with ZR and LEA"),
    };
    PlotOptions optComparison;
    optComparison.title = "Total comparison";
    optComparison.canvasName = "cComparison";
    optComparison.saveAs = "./Figures/compare_comparison.png";
    DrawComparison(calc_comparison, optComparison);

    // --- Total comparison  with ZR LEA ---
    std::vector<Calc> calc_comparison_ZR_LEA = {
        Calc("../../FranckCalculations/DWUCK4_personal_calculations/good_calculations/v1/xs_dw4.n", "DWUCK4",
             Format::DWUCK4, AutoColor(0), 1)
            .Tag("both SO - my version 3 MeV"),
        Calc("../../FranckCalculations/DWUCK4_examples_Franck/3MeV_Ex_pwave_noOscilation/xs_dw4.n", "DWUCK4",
             Format::DWUCK4, AutoColor(1), 1)
            .Tag("no SO - Franck's version 3 MeV"),
        Calc("./7Lidp8Li_prior/fort.202", "FRESCO", Format::FRESCO, AutoColor(3), 1)
            .Tag("ADWA prior - Mario's version near 3 MeV"),
        Calc("./7Lidp8Li_DWBA_prior/fort.202", "FRESCO", Format::FRESCO, AutoColor(5), 1)
            .Tag("DWBA prior - Mario's version near 3 MeV"),
        Calc("./personal_calculations/DWBA/prior/ALL/fort.202", "FRESCO", Format::FRESCO, AutoColor(8), 1)
            .Tag("DWBA prior - All changes with ZR and LEA"),
        Calc("./personal_calculations/ADWA/prior/ALL/fort.202", "FRESCO", Format::FRESCO, AutoColor(9), 1)
            .Tag("ADWA prior - All changes with ZR and LEA"),
    };
    PlotOptions optComparison_ZR_LEA;
    optComparison_ZR_LEA.title = "Total comparison ZR LEA";
    optComparison_ZR_LEA.canvasName = "cComparison ZR LEA";
    optComparison_ZR_LEA.saveAs = "./Figures/compare_comparison_ZR_LEA.png";
    DrawComparison(calc_comparison_ZR_LEA, optComparison_ZR_LEA);
}