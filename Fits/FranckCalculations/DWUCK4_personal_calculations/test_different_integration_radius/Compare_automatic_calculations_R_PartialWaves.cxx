// =============================================================================
//  Compare_automatic_calculations_R_PartialWaves.C
//
//  Dibuja las distribuciones angulares de DWUCK4 (xs_dw4.n) producidas por
//  MakeDwuck_R_PartialWaves.C con la estructura:
//        <step>/R<R>/L<L>/xs_dw4.n        p.ej.  S01/R40/L30/xs_dw4.n
//  La etiqueta del step (DRF) sale de la tarjeta 4 (linea 6) de run.sh.
//
//  Figuras (en ./Figures), a un Lmax fijo "lPlot":
//     compare_<step>_L<L>.png    un canvas por step: todos los R (color = R)
//     ratio_<step>_L<L>.png      por step: sigma(R)/sigma(R mayor del step)
//     compare_total_L<L>.png     TOTAL: todos los steps y R (color = R, linea = step)
//     ratio_total_L<L>.png       TOTAL vs referencia global (mayor R, menor step)
//     lconv_<step>_R<R>.png      (opcional) convergencia en L de cada (step,R)
//
//  Uso:
//     root -l 'Compare_automatic_calculations_R_PartialWaves.C+'
//     root -l 'Compare_automatic_calculations_R_PartialWaves.C+({"S01","S02"},{40,80,120},{30,60},60)'
//     root -l 'Compare_automatic_calculations_R_PartialWaves.C+({"S01","S02","S03"},{},{},0,true,"E_{x} = 3 MeV",0.1,true)'
//
//  Si falta algun archivo se avisa y se omite esa curva.
// =============================================================================

#include "Rtypes.h"

#include <TAxis.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMultiGraph.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace CompDW
{

// ----------------------------- CONFIGURACION ---------------------------------
const char* const kXsFile = "xs_dw4.n";
const char* const kRunFile = "run.sh";
const char* const kFigDir = "./Figures";

const char* const kXTitle = "#theta_{cm} (deg)";
const char* const kYTitle = "d#sigma/d#Omega (mb/sr)";

inline std::vector<int> Range(int a, int b, int s)
{
    std::vector<int> v;
    for(int x = a; x <= b; x += s)
        v.push_back(x);
    return v;
}

// ------------------------------ utilidades -----------------------------------
struct Curve
{
    TString label;
    Color_t color;
    Style_t style;
    std::vector<double> x, y;
};

inline Color_t ColorIdx(size_t i)
{
    static const Color_t pal[] = {kBlack,       kRed + 1,    kBlue + 1,  kGreen + 2,  kMagenta + 1, kOrange + 7,
                                  kCyan + 2,    kViolet + 1, kYellow + 2, kPink + 7,  kTeal + 3,    kSpring + 4,
                                  kAzure + 1,   kGray + 2};
    return pal[i % (sizeof(pal) / sizeof(pal[0]))];
}

inline Style_t StyleIdx(size_t i)
{
    static const Style_t sty[] = {1, 2, 7, 3, 9, 5};
    return sty[i % (sizeof(sty) / sizeof(sty[0]))];
}

// Lee "angulo xs" (acepta exponentes Fortran 1.0D-03). Lineas no numericas se ignoran.
inline bool ReadXs(const TString& path, std::vector<double>& x, std::vector<double>& y)
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
        x.push_back(a);
        y.push_back(b);
    }
    return !x.empty();
}

inline Curve Load(const TString& path, const TString& label, Color_t col, Style_t sty, bool warn = true)
{
    Curve c;
    c.label = label;
    c.color = col;
    c.style = sty;
    if(!ReadXs(path, c.x, c.y))
    {
        if(warn)
            std::cerr << "[AVISO] sin datos / no existe: " << path << " (se omite)" << std::endl;
        c.x.clear();
        c.y.clear();
    }
    return c;
}

inline TString XsPath(const std::string& step, int R, int L)
{
    return Form("%s/R%d/L%d/%s", step.c_str(), R, L, kXsFile);
}

// DRF (cols 1-8) y RMAX (cols 17-24) de la tarjeta 4 = linea 6 de un run.sh
inline bool ReadCard4(const TString& runsh, double& drf, double& rmax)
{
    std::ifstream in(runsh.Data());
    if(!in)
        return false;
    std::string line;
    for(int n = 1; std::getline(in, line); ++n)
    {
        if(n != 6)
            continue;
        if(line.size() < 24)
            line.resize(24, ' ');
        drf = std::atof(line.substr(0, 8).c_str());
        rmax = std::atof(line.substr(16, 8).c_str());
        return drf > 0. && rmax > 0.;
    }
    return false;
}

// DRF del step leido de cualquier run.sh existente; -1 si no hay ninguno
inline double StepDRF(const std::string& step, const std::vector<int>& rList, const std::vector<int>& lList)
{
    for(size_t ir = 0; ir < rList.size(); ++ir)
        for(size_t j = 0; j < lList.size(); ++j)
        {
            double drf = 0., rmax = 0.;
            if(ReadCard4(Form("%s/R%d/L%d/%s", step.c_str(), rList[ir], lList[j], kRunFile), drf, rmax))
                return drf;
        }
    return -1.;
}

// sigma(i)/sigma(iref); solo donde sigma(ref) > 0. La curva de referencia se omite.
inline std::vector<Curve> Ratios(const std::vector<Curve>& cv, size_t iref)
{
    std::vector<Curve> out;
    if(iref >= cv.size() || cv[iref].x.empty())
    {
        std::cerr << "[AVISO] no hay curva de referencia: no se dibujan cocientes" << std::endl;
        return out;
    }
    const Curve& ref = cv[iref];
    for(size_t i = 0; i < cv.size(); ++i)
    {
        if(i == iref || cv[i].x.empty())
            continue;
        if(cv[i].x.size() != ref.x.size() || std::fabs(cv[i].x.front() - ref.x.front()) > 1e-6 ||
           std::fabs(cv[i].x.back() - ref.x.back()) > 1e-6)
        {
            std::cerr << "[AVISO] " << cv[i].label << ": distinta malla angular que la referencia (se omite)"
                      << std::endl;
            continue;
        }
        Curve r;
        r.label = cv[i].label;
        r.color = cv[i].color;
        r.style = cv[i].style;
        for(size_t k = 0; k < ref.x.size(); ++k)
            if(ref.y[k] > 0.)
            {
                r.x.push_back(ref.x[k]);
                r.y.push_back(cv[i].y[k] / ref.y[k]);
            }
        out.push_back(r);
    }
    return out;
}

// ------------------------- canvas con una lista de curvas ---------------------
inline TCanvas* DrawCurves(const std::vector<Curve>& cv, const TString& title, const TString& name,
                           const TString& saveAs, bool logY, const TString& yTitle, double ymin = 0., double ymax = 0.,
                           bool lineAtOne = false)
{
    gStyle->SetOptStat(0);

    auto* c = new TCanvas(name, name, 900, 650);
    c->SetLeftMargin(0.13);
    c->SetGrid();
    if(logY)
        c->SetLogy();

    auto* mg = new TMultiGraph(Form("mg_%s", name.Data()), "");

    std::vector<TGraph*> gs;
    for(size_t i = 0; i < cv.size(); ++i)
    {
        std::vector<double> x, y;
        for(size_t k = 0; k < cv[i].x.size(); ++k)
            if(!logY || cv[i].y[k] > 0.)
            {
                x.push_back(cv[i].x[k]);
                y.push_back(cv[i].y[k]);
            }
        if(x.empty())
            continue;
        auto* g = new TGraph(static_cast<int>(x.size()), x.data(), y.data());
        g->SetName(Form("g_%s_%d", name.Data(), static_cast<int>(i)));
        g->SetTitle(cv[i].label);
        g->SetLineColor(cv[i].color);
        g->SetLineStyle(cv[i].style);
        g->SetLineWidth(2);
        g->SetMarkerStyle(0);
        g->SetFillStyle(0);
        mg->Add(g, "L");
        gs.push_back(g);
    }

    if(gs.empty())
    {
        std::cerr << "[ERROR] " << name << ": ningun archivo pudo leerse, no hay nada que dibujar." << std::endl;
        return c;
    }

    const int n = static_cast<int>(gs.size());
    const int ncol = (n > 16) ? 3 : (n > 8) ? 2 : 1;
    const int rows = (n + ncol - 1) / ncol;
    const double wcol = (ncol == 3) ? 0.22 : 0.27;
    auto* leg = new TLegend(0.90 - wcol * ncol, 0.88 - 0.04 * rows, 0.90, 0.88);
    leg->SetNColumns(ncol);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    if(n > 8)
        leg->SetTextSize(0.027);
    for(size_t i = 0; i < gs.size(); ++i)
        leg->AddEntry(gs[i], gs[i]->GetTitle(), "L");

    mg->SetTitle(Form("%s;%s;%s", title.Data(), kXTitle, yTitle.Data()));
    if(ymax > ymin)
    {
        mg->SetMinimum(ymin);
        mg->SetMaximum(ymax);
    }
    mg->Draw("A");
    mg->GetXaxis()->SetTitleOffset(1.1);
    mg->GetYaxis()->SetTitleOffset(1.3);
    if(lineAtOne)
    {
        auto* l = new TLine(mg->GetXaxis()->GetXmin(), 1., mg->GetXaxis()->GetXmax(), 1.);
        l->SetLineColor(kGray + 1);
        l->SetLineStyle(1);
        l->Draw();
    }
    leg->Draw();

    c->Modified();
    c->Update();
    if(saveAs.Length() > 0)
        c->SaveAs(saveAs);
    return c;
}

} // namespace CompDW

// -----------------------------------------------------------------------------
//  Punto de entrada
//     stepDirs   carpetas de step (las de MakeDwuck_R_PartialWaves.C)
//     rList      RMAX a dibujar (vacio -> 40,50,...,120)
//     lList      Lmax disponibles (vacio -> 30,35,...,60); se usan para leer run.sh y para lconv
//     lPlot      Lmax con el que se dibujan los canvas por step y el total (0 -> mayor L con datos)
//     logY       escala log en las distribuciones (los cocientes son siempre lineales)
//     exLabel    texto del titulo
//     ratioRange los cocientes se muestran en 1 +- ratioRange
//     drawLconv  si true, ademas un canvas de convergencia en L por cada (step,R)
// -----------------------------------------------------------------------------
void Compare_automatic_calculations_R_PartialWaves(const std::vector<std::string>& stepDirs = {"S01", "S02", "S03"},
                                                   const std::vector<int>& rIn = std::vector<int>(),
                                                   const std::vector<int>& lIn = std::vector<int>(), int lPlot = 30,
                                                   bool logY = true, const std::string& exLabel = "E_{x} = 3 MeV",
                                                   double ratioRange = 0.2, bool drawLconv = false)
{
    using namespace CompDW;

    const std::vector<int> rList = rIn.empty() ? Range(40, 120, 10) : rIn;
    const std::vector<int> lList = lIn.empty() ? Range(30, 60, 5) : lIn;

    if(stepDirs.empty())
    {
        std::cerr << "[ERROR] lista de steps vacia" << std::endl;
        return;
    }
    gSystem->mkdir(kFigDir, true);

    const size_t ns = stepDirs.size(), nr = rList.size();

    // ---- Lmax de trabajo: el pedido, o el mayor L para el que exista algun fichero ----
    if(lPlot <= 0)
    {
        for(size_t j = 0; j < lList.size(); ++j)
        {
            bool any = false;
            for(size_t is = 0; is < ns && !any; ++is)
                for(size_t ir = 0; ir < nr && !any; ++ir)
                {
                    std::ifstream f(XsPath(stepDirs[is], rList[ir], lList[j]).Data());
                    any = bool(f);
                }
            if(any && lList[j] > lPlot)
                lPlot = lList[j];
        }
        if(lPlot <= 0)
        {
            std::cerr << "[ERROR] no se encontro ningun xs_dw4.n" << std::endl;
            return;
        }
    }
    std::cout << "Lmax usado en los canvas por step y en el total: " << lPlot << std::endl;

    // ---- curvas a lPlot: curves[is][ir] ----
    std::vector<std::vector<Curve>> curves(ns);
    std::vector<TString> stepLab(ns);
    for(size_t is = 0; is < ns; ++is)
    {
        const double drf = StepDRF(stepDirs[is], rList, lList);
        stepLab[is] = drf > 0. ? TString(Form("%s (#Deltar=%g fm)", stepDirs[is].c_str(), drf))
                               : TString(stepDirs[is].c_str());
        for(size_t ir = 0; ir < nr; ++ir)
            curves[is].push_back(
                Load(XsPath(stepDirs[is], rList[ir], lPlot), Form("R_{max}=%d fm", rList[ir]), ColorIdx(ir), 1));
    }

    // indice del mayor R con datos dentro de un step; -1 si no hay
    auto refR = [&](size_t is) {
        int jr = -1;
        for(size_t ir = 0; ir < nr; ++ir)
            if(!curves[is][ir].x.empty() && (jr < 0 || rList[ir] > rList[static_cast<size_t>(jr)]))
                jr = static_cast<int>(ir);
        return jr;
    };

    // ---- 1) un canvas por step (y sus cocientes respecto al mayor R) ----
    for(size_t is = 0; is < ns; ++is)
    {
        const TString title = TString(exLabel.c_str()) + ", " + stepLab[is] + Form(", L_{max}=%d", lPlot);
        const TString tag = Form("%s_L%d", stepDirs[is].c_str(), lPlot);

        DrawCurves(curves[is], title, TString("cStep_") + tag, TString(kFigDir) + "/compare_" + tag + ".png", logY,
                   kYTitle);

        const int jr = refR(is);
        if(jr < 0)
            continue;
        std::vector<Curve> rat = Ratios(curves[is], static_cast<size_t>(jr));
        DrawCurves(rat, title + Form(" (ref. R_{max}=%d)", rList[static_cast<size_t>(jr)]),
                   TString("cStepRatio_") + tag, TString(kFigDir) + "/ratio_" + tag + ".png", false,
                   Form("d#sigma(R_{max}) / d#sigma(R_{max}=%d)", rList[static_cast<size_t>(jr)]), 1. - ratioRange,
                   1. + ratioRange, true);
    }

    // ---- 2) TOTAL: todos los steps y todos los R (color = R, tipo de linea = step) ----
    {
        std::vector<Curve> all;
        int iref = -1; // referencia global: mayor R, y a igualdad menor DRF (= primer step de la lista con ese R)
        for(size_t ir = 0; ir < nr; ++ir)
            for(size_t is = 0; is < ns; ++is)
            {
                Curve c = curves[is][ir];
                c.label = Form("%s, R_{max}=%d", stepDirs[is].c_str(), rList[ir]);
                c.color = ColorIdx(ir);
                c.style = StyleIdx(is);
                all.push_back(c);
            }
        DrawCurves(all, TString(exLabel.c_str()) + Form(", L_{max}=%d", lPlot), "cTotal",
                   TString(kFigDir) + Form("/compare_total_L%d.png", lPlot), logY, kYTitle);

        // referencia: mayor R con datos; entre los steps, el de menor DRF
        double bestDrf = 1e9;
        int bestR = -1;
        for(size_t is = 0; is < ns; ++is)
        {
            const int jr = refR(is);
            if(jr < 0)
                continue;
            const int R = rList[static_cast<size_t>(jr)];
            const double drf = StepDRF(stepDirs[is], rList, lList);
            const double d = drf > 0. ? drf : 1e8;
            if(R > bestR || (R == bestR && d < bestDrf))
            {
                bestR = R;
                bestDrf = d;
                // posicion en 'all' (R externo, step interno)
                iref = static_cast<int>(static_cast<size_t>(jr) * ns + is);
            }
        }
        if(iref >= 0)
        {
            const std::string refName = all[static_cast<size_t>(iref)].label.Data();
            std::vector<Curve> rat = Ratios(all, static_cast<size_t>(iref));
            DrawCurves(rat, TString(exLabel.c_str()) + Form(", L_{max}=%d", lPlot) + " (ref. " + refName.c_str() + ")",
                       "cTotalRatio", TString(kFigDir) + Form("/ratio_total_L%d.png", lPlot), false,
                       TString("d#sigma / d#sigma(") + refName.c_str() + ")", 1. - ratioRange, 1. + ratioRange,
                       true);
        }
    }

    // ---- 3) opcional: convergencia en L de cada (step,R) ----
    if(drawLconv)
    {
        for(size_t is = 0; is < ns; ++is)
            for(size_t ir = 0; ir < nr; ++ir)
            {
                std::vector<Curve> cl;
                int jref = -1;
                for(size_t j = 0; j < lList.size(); ++j)
                {
                    cl.push_back(Load(XsPath(stepDirs[is], rList[ir], lList[j]), Form("L_{max}=%d", lList[j]),
                                      ColorIdx(j), 1, false));
                    if(!cl.back().x.empty() && (jref < 0 || lList[j] > lList[static_cast<size_t>(jref)]))
                        jref = static_cast<int>(j);
                }
                if(jref < 0)
                    continue;
                const int lRef = lList[static_cast<size_t>(jref)];
                const TString title = TString(exLabel.c_str()) + ", " + stepLab[is] + Form(", R_{max}=%d fm", rList[ir]);
                const TString tag = Form("%s_R%d", stepDirs[is].c_str(), rList[ir]);
                std::vector<Curve> rat = Ratios(cl, static_cast<size_t>(jref));
                DrawCurves(rat, title + Form(" (ref. L_{max}=%d)", lRef), TString("cLconv_") + tag,
                           TString(kFigDir) + "/lconv_" + tag + ".png", false,
                           Form("d#sigma(L_{max}) / d#sigma(L_{max}=%d)", lRef), 1. - ratioRange, 1. + ratioRange,
                           true);
            }
    }
}