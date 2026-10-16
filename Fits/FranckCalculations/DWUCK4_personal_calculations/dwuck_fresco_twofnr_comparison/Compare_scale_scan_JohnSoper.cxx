// =============================================================================
//  Compare_scale_scan.C
//
//  Macro de ROOT: recorre las carpetas DifferentScaleFactor/VSOR<x>_VSOI<y>/
//  generadas por el script de escalado, compara cada xs_dw4.n con la seccion
//  eficaz de FRESCO, ordena las combinaciones por diferencia y dibuja FRESCO
//  junto con las N mejores.
//
//  Metrica (por combinacion):
//    d_i   = ln( sigma_DWUCK(theta_i) / sigma_FRESCO(theta_i) )
//    score = sqrt( mean( d_i^2 ) )          (RMS del log del cociente)
//  - DWUCK4 se interpola (lineal en log) en los angulos de FRESCO.
//  - score ~ 0.10  ->  ~10 % de diferencia tipica entre ambas curvas.
//  - shapeOnly = true: se resta el cociente medio (compara solo la FORMA,
//    ignorando la normalizacion global).
//
//  Filtros de combinaciones (antes de leer ningun fichero):
//    sameValue = true  -> solo parejas con |VSOR| == |VSOI|
//    signVSOR / signVSOI = Any | Positive | Negative
//    (ejemplo: VSOR positivo con VSOI negativo -> Positive / Negative)
//
//  Uso:
//    root -l Compare_scale_scan.C
//    root -l 'Compare_scale_scan.C+'      // compilada
// =============================================================================

#include "Rtypes.h"

#include <TAxis.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
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

struct Calc
{
    TString path;
    TString label;
    Format format;
    Color_t color; // -1 => automatico
    Style_t lineStyle;
    Width_t lineWidth;
    TString tag;

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

    Calc& Tag(const TString& t)
    {
        tag = t;
        return *this;
    }

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
            continue;
        if(useEnd && IsEndMarker(line))
            break;
        if(Trim(line).empty())
            continue;

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
//  Construccion del TGraph
// -----------------------------------------------------------------------------

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

inline Color_t AutoColor(size_t i)
{
    static const Color_t palette[] = {kBlack,      kRed + 1,  kBlue + 1,   kGreen + 2, kMagenta + 1,
                                      kOrange + 7, kCyan + 2, kViolet + 1, kGray + 2,  kYellow + 2};
    return palette[i % (sizeof(palette) / sizeof(palette[0]))];
}

inline TCanvas* DrawComparison(const std::vector<Calc>& calcs, const PlotOptions& opt = PlotOptions())
{
    gStyle->SetOptStat(0);

    auto* c = new TCanvas(opt.canvasName, opt.canvasName, opt.width, opt.height);
    c->SetLeftMargin(0.13);
    c->SetGrid();
    if(opt.logY)
        c->SetLogy();

    auto* mg = new TMultiGraph(Form("mg_%s", opt.canvasName.Data()), "");
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

// =============================================================================
//  NUEVO: barrido de factores de escala y ranking frente a FRESCO
// =============================================================================

// Filtro de signo para una de las dos variables
enum class SignSel
{
    Any,      // sin filtro
    Positive, // solo valores > 0
    Negative  // solo valores < 0
};

struct ScanOptions
{
    TString scanDir = "DifferentScaleFactor"; // carpeta que contiene VSOR*_VSOI*
    TString xsFile = "xs_dw4.n";              // fichero dentro de cada subcarpeta
    std::vector<TString> vsorList;            // TEXTO exacto usado en el nombre de carpeta
    std::vector<TString> vsoiList;
    double thetaMin = 0.0; // rango angular usado en la comparacion numerica
    double thetaMax = 180.0;
    int nBest = 5;                                     // cuantas curvas dibujar
    bool shapeOnly = false;                            // true: ignora la normalizacion global
    bool drawMap = true;                               // mapa 2D del score (VSOR vs VSOI)
    TString rankingFile = "ranking_scale_factors.txt"; // "" para no escribir

    // --- Filtros de combinaciones ---
    bool sameValue = false;          // true: solo parejas con |VSOR| == |VSOI|
    SignSel signVSOR = SignSel::Any; // signo permitido para VSOR
    SignSel signVSOI = SignSel::Any; // signo permitido para VSOI
};

struct ScanResult
{
    TString vsor, vsoi, path;
    double score = 0; // RMS del ln(DWUCK/FRESCO)
    double bias = 0;  // media del ln(DWUCK/FRESCO)
    int n = 0;        // puntos usados
};

inline bool SignOk(double v, SignSel s)
{
    switch(s)
    {
    case SignSel::Positive: return v > 0;
    case SignSel::Negative: return v < 0;
    default: break;
    }
    return true;
}

inline const char* SignName(SignSel s)
{
    return s == SignSel::Positive ? "+" : (s == SignSel::Negative ? "-" : "any");
}

// Decide si una pareja (VSOR, VSOI) pasa los filtros. Los valores se leen del TEXTO de las listas.
inline bool PairSelected(const TString& vsor, const TString& vsoi, const ScanOptions& so)
{
    const double a = vsor.Atof();
    const double b = vsoi.Atof();
    if(!SignOk(a, so.signVSOR) || !SignOk(b, so.signVSOI))
        return false;
    if(so.sameValue)
    {
        const double ma = std::fabs(a), mb = std::fabs(b);
        if(std::fabs(ma - mb) > 1e-9 * std::max(1.0, std::max(ma, mb)))
            return false;
    }
    return true;
}

// Interpolacion lineal en log(y). x debe estar ordenado de menor a mayor.
inline bool InterpLogY(const std::vector<double>& x, const std::vector<double>& y, double xq, double& yq)
{
    if(x.size() < 2 || xq < x.front() || xq > x.back())
        return false;
    auto it = std::upper_bound(x.begin(), x.end(), xq);
    size_t i1 = static_cast<size_t>(it - x.begin());
    if(i1 >= x.size())
        i1 = x.size() - 1;
    if(i1 == 0)
        return false;
    const size_t i0 = i1 - 1;
    if(y[i0] <= 0 || y[i1] <= 0)
        return false;
    const double t = (x[i1] > x[i0]) ? (xq - x[i0]) / (x[i1] - x[i0]) : 0.0;
    yq = std::exp(std::log(y[i0]) * (1.0 - t) + std::log(y[i1]) * t);
    return true;
}

// Devuelve false si hay menos de 3 puntos comparables.
inline bool CompareCurves(const std::vector<double>& frx, const std::vector<double>& fry,
                          const std::vector<double>& dwx, const std::vector<double>& dwy, const ScanOptions& so,
                          ScanResult& res)
{
    std::vector<double> d;
    for(size_t i = 0; i < frx.size(); ++i)
    {
        if(frx[i] < so.thetaMin || frx[i] > so.thetaMax || fry[i] <= 0)
            continue;
        double yq;
        if(!InterpLogY(dwx, dwy, frx[i], yq))
            continue;
        d.push_back(std::log(yq / fry[i]));
    }
    if(d.size() < 3)
        return false;

    double mean = 0;
    for(double v : d)
        mean += v;
    mean /= d.size();

    double ss = 0;
    for(double v : d)
    {
        const double w = so.shapeOnly ? (v - mean) : v;
        ss += w * w;
    }
    res.score = std::sqrt(ss / d.size());
    res.bias = mean;
    res.n = static_cast<int>(d.size());
    return true;
}

// Recorre todas las combinaciones, las ordena y dibuja FRESCO + las N mejores.
inline std::vector<ScanResult>
ScanScaleFactors(const Calc& fresco, const ScanOptions& so, const PlotOptions& po = PlotOptions())
{
    std::vector<ScanResult> results;

    std::vector<double> frx, fry;
    if(!ReadData(fresco, frx, fry))
        return results;

    int nMissing = 0, nUnreadable = 0, nNoOverlap = 0, nSelected = 0;
    for(const auto& vsor : so.vsorList)
        for(const auto& vsoi : so.vsoiList)
        {
            if(!PairSelected(vsor, vsoi, so))
                continue;
            ++nSelected;

            ScanResult r;
            r.vsor = vsor;
            r.vsoi = vsoi;
            r.path = Form("%s/VSOR%s_VSOI%s/%s", so.scanDir.Data(), vsor.Data(), vsoi.Data(), so.xsFile.Data());

            if(gSystem->AccessPathName(r.path)) // true si NO existe
            {
                ++nMissing;
                continue;
            }
            std::vector<double> dwx, dwy;
            Calc dw(r.path, "", Format::DWUCK4);
            if(!ReadData(dw, dwx, dwy))
            {
                ++nUnreadable;
                continue;
            }
            if(!CompareCurves(frx, fry, dwx, dwy, so, r))
            {
                ++nNoOverlap;
                continue;
            }
            results.push_back(r);
        }

    const int nAll = static_cast<int>(so.vsorList.size() * so.vsoiList.size());
    std::cout << "\n[SCAN] Filtro: sameValue=" << (so.sameValue ? "si" : "no")
              << ", signo VSOR=" << SignName(so.signVSOR) << ", signo VSOI=" << SignName(so.signVSOI) << "  ->  "
              << nSelected << " de " << nAll << " combinaciones seleccionadas" << std::endl;
    std::cout << "[SCAN] " << results.size() << " / " << nSelected << " evaluadas"
              << " (faltan: " << nMissing << ", ilegibles: " << nUnreadable << ", sin solape angular: " << nNoOverlap
              << ")" << std::endl;
    if(results.empty())
    {
        if(nSelected == 0)
            std::cerr << "[ERROR] Ninguna combinacion pasa el filtro. Con sameValue=true ambas listas deben "
                         "contener el mismo valor absoluto, y el signo elegido debe aparecer en cada lista."
                      << std::endl;
        else
            std::cerr << "[ERROR] No hay resultados. Revisa scanDir y el texto de las listas (deben coincidir con los "
                         "nombres de carpeta)."
                      << std::endl;
        return results;
    }

    std::sort(results.begin(), results.end(),
              [](const ScanResult& a, const ScanResult& b) { return a.score < b.score; });

    // --- Tabla por pantalla (mejores) -------------------------------------------------
    const int nShow = std::min<int>(std::max(so.nBest, 10), static_cast<int>(results.size()));
    std::cout << "\n  rank   VSOR       VSOI       score    norm.factor   n\n";
    std::cout << "  ---------------------------------------------------------\n";
    for(int i = 0; i < nShow; ++i)
        std::printf("  %3d   %-9s  %-9s  %7.4f  %10.4f   %d\n", i + 1, results[i].vsor.Data(), results[i].vsoi.Data(),
                    results[i].score, std::exp(-results[i].bias), results[i].n);
    std::cout << "\n  norm.factor = factor por el que habria que multiplicar DWUCK4 para igualar la\n"
                 "  normalizacion media de FRESCO (1.0 = ya coinciden en magnitud).\n"
              << std::endl;

    // --- Fichero con TODAS las combinaciones ---------------------------------------
    if(so.rankingFile.Length() > 0)
    {
        std::ofstream out(so.rankingFile.Data());
        out << "# rank  VSOR  VSOI  score  norm_factor  npts   (shapeOnly=" << so.shapeOnly << ", theta in ["
            << so.thetaMin << "," << so.thetaMax << "], sameValue=" << so.sameValue
            << ", signVSOR=" << SignName(so.signVSOR) << ", signVSOI=" << SignName(so.signVSOI) << ")\n";
        for(size_t i = 0; i < results.size(); ++i)
            out << i + 1 << "  " << results[i].vsor << "  " << results[i].vsoi << "  " << results[i].score << "  "
                << std::exp(-results[i].bias) << "  " << results[i].n << "\n";
        std::cout << "[SCAN] Ranking completo guardado en " << so.rankingFile << std::endl;
    }

    // --- Figura: FRESCO + N mejores ---------------------------------------------------
    std::vector<Calc> toPlot;
    Calc ref = fresco;
    if(ref.color < 0)
        ref.color = kBlack;
    toPlot.push_back(ref);

    const int nPlot = std::min<int>(so.nBest, static_cast<int>(results.size()));
    for(int i = 0; i < nPlot; ++i)
    {
        Calc c(results[i].path, "DWUCK4", Format::DWUCK4, AutoColor(i + 1), 1, 2);
        c.Tag(Form("VSOR=%s, VSOI=%s, rms=%.3f", results[i].vsor.Data(), results[i].vsoi.Data(), results[i].score));
        toPlot.push_back(c);
    }
    DrawComparison(toPlot, po);

    // --- Mapa 2D del score ------------------------------------------------------------
    if(so.drawMap)
    {
        const int nR = static_cast<int>(so.vsorList.size());
        const int nI = static_cast<int>(so.vsoiList.size());
        delete gROOT->FindObject("hScoreMap");
        auto* h = new TH2D("hScoreMap", Form("score (RMS ln DWUCK/FRESCO);VSOR;VSOI"), nR, 0, nR, nI, 0, nI);
        h->SetStats(0);
        for(int i = 0; i < nR; ++i)
            h->GetXaxis()->SetBinLabel(i + 1, so.vsorList[i]);
        for(int j = 0; j < nI; ++j)
            h->GetYaxis()->SetBinLabel(j + 1, so.vsoiList[j]);
        for(const auto& r : results)
        {
            int ix = -1, iy = -1;
            for(int i = 0; i < nR; ++i)
                if(so.vsorList[i] == r.vsor)
                    ix = i;
            for(int j = 0; j < nI; ++j)
                if(so.vsoiList[j] == r.vsoi)
                    iy = j;
            if(ix >= 0 && iy >= 0)
                h->SetBinContent(ix + 1, iy + 1, r.score);
        }
        h->SetMinimum(results.front().score * 0.999); // las combinaciones sin dato (0) quedan en blanco
        gStyle->SetPaintTextFormat(".2f");
        h->SetMarkerSize(0.9);
        auto* cm = new TCanvas("cScoreMap", "cScoreMap", 900, 700);
        cm->SetRightMargin(0.15);
        h->Draw("COLZ TEXT");
        cm->Modified();
        cm->Update();
    }

    return results;
}

} // namespace CompCalc

// -----------------------------------------------------------------------------
//  Punto de entrada: EDITA AQUI
// -----------------------------------------------------------------------------
void Compare_scale_scan_JohnSoper()
{
    using namespace CompCalc;

    // FRESCO (referencia)
    Calc fresco("../../FRESCO_personal_calculations/calcs_from_personal_twofnr_calcs/"
                "gs_ADWA_KD_twofnr_JohnSoper_LEA(zerorange)_noSO_exit/fort.202",
                "FRESCO", Format::FRESCO, kBlack, 1, 3);
    fresco.Tag("ADKD-ZR-LEA, noSO exit");

    ScanOptions so;
    // Carpeta donde el script de escalado creo DifferentScaleFactor/
    so.scanDir = "./gs_ADWA_KD_twofnr_JohnSoper_LEA(zerorange)_noSO_exit_diferentScale_soEntrance/"
                 "DifferentScaleFactor";

    // IMPORTANTE: el texto debe ser identico al usado en el script bash
    // (carpeta VSOR-452_VSOI4520  ->  "-452" y "4520")
    // Lista ampliada (12x12 = 144 combinaciones):
    so.vsorList = {"-1", "1", "45.2", "-45.2", "22.6", "-22.6", "11.3", "-11.3", "90.4", "-90.4", "180.8", "-180.8"};
    so.vsoiList = {"-1", "1", "45.2", "-45.2", "22.6", "-22.6", "11.3", "-11.3", "90.4", "-90.4", "180.8", "-180.8"};
    // Lista actual (10x10 = 100 combinaciones):
    // so.vsorList = {"-1", "1", "45.2", "-45.2", "22.6", "-22.6", "11.3", "-11.3", "90.4", "-90.4"};
    // so.vsoiList = {"-1", "1", "45.2", "-45.2", "22.6", "-22.6", "11.3", "-11.3", "90.4", "-90.4"};

    // --- Filtros de combinaciones ---
    // sameValue: true -> solo parejas con el mismo valor absoluto (|VSOR| == |VSOI|).
    //            Ambas listas deben contener ese valor (p.ej. "452" en las dos).
    so.sameValue = true;
    // Signo permitido: SignSel::Any (sin filtro), SignSel::Positive o SignSel::Negative.
    //   VSOR+ con VSOI+ : Positive / Positive      VSOR+ con VSOI- : Positive / Negative
    //   VSOR- con VSOI+ : Negative / Positive      VSOR- con VSOI- : Negative / Negative
    so.signVSOR = SignSel::Any;
    so.signVSOI = SignSel::Any;

    so.thetaMin = 0.0;  // rango angular de la comparacion numerica
    so.thetaMax = 80.0; // ajusta a la region que te interese
    so.nBest = 5;
    so.shapeOnly = false; // true: comparar solo la forma, ignorando la normalizacion

    PlotOptions po;
    po.title = "g.s. no SO exit: best scale factors vs FRESCO";
    po.canvasName = "cBestScale";
    po.saveAs = "./Figures/compare_gs_bestscale.png";

    ScanScaleFactors(fresco, so, po);
}