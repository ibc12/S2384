#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TStyle.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// --- Estructura para guardar parámetros ---
struct CalData {
    std::vector<double> a_vals;
    std::vector<double> b_vals;
};

// --- Leer archivo ---
CalData ReadCalibrationFileRaw(const char* filename)
{
    std::ifstream in(filename);
    if (!in.is_open()) {
        std::cerr << "Error abriendo archivo " << filename << std::endl;
        return {};
    }

    std::vector<double> a_vals, b_vals;
    std::string name;
    double a, b;

    while (in >> name >> a >> b) {
        if (name.find("_E") != std::string::npos) {
            a_vals.push_back(a);
            b_vals.push_back(b);
        }
    }
    in.close();
    return {a_vals, b_vals};
}

// --- Graficos de parametros ---
TGraph* MakeGraph(const std::vector<double>& vals, int color, int markerStyle)
{
    int n = vals.size();
    std::vector<double> idx(n);
    for (int i = 0; i < n; i++) idx[i] = i;

    TGraph* g = new TGraph(n, idx.data(), vals.data());
    g->SetMarkerStyle(markerStyle);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    return g;
}

// --- Grafico de conversion canal 3000 ---
TGraph* MakeGraphConv(const CalData& data, int color, int markerStyle, double channel = 3000)
{
    int n = data.a_vals.size();
    std::vector<double> idx(n), conv(n);
    for (int i = 0; i < n; i++) {
        idx[i] = i;
        conv[i] = data.a_vals[i] + data.b_vals[i]*channel;
    }
    TGraph* g = new TGraph(n, idx.data(), conv.data());
    g->SetMarkerStyle(markerStyle);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    return g;
}

// --- Dibujar canvas con 3 pads: a, b, E(3000) ---
void DrawComparison(const CalData& pre, const CalData& post,
                    const char* title, const char* tagPre, const char* tagPost)
{
    TCanvas* c = new TCanvas(title, title, 1000, 800);
    c->Divide(1,3);

    // Pad 1: a
    c->cd(1);
    auto gA_pre  = MakeGraph(pre.a_vals, kBlue, 21);
    auto gA_post = MakeGraph(post.a_vals, kRed, 20);
    gA_post->SetTitle(TString::Format("%s - Parametro a;Indice;a", title));
    gA_post->Draw("APL");
    gA_pre->Draw("PL SAME");
    auto leg1 = new TLegend(0.7,0.75,0.9,0.9);
    leg1->AddEntry(gA_post, tagPost, "lp");
    leg1->AddEntry(gA_pre, tagPre, "lp");
    leg1->Draw();

    // Pad 2: b
    c->cd(2);
    auto gB_pre  = MakeGraph(pre.b_vals, kBlue, 21);
    auto gB_post = MakeGraph(post.b_vals, kRed, 20);
    gB_post->SetTitle(TString::Format("%s - Parametro b;Indice;b", title));
    gB_post->Draw("APL");
    gB_pre->Draw("PL SAME");

    // Pad 3: E(3000)
    c->cd(3);
    auto gE_pre  = MakeGraphConv(pre, kBlue, 21);
    auto gE_post = MakeGraphConv(post, kRed, 20);
    gE_post->SetTitle(TString::Format("%s - E(canal=3000);Indice;E [keV]", title));
    gE_post->Draw("APL");
    gE_pre->Draw("PL SAME");
}

// --- Dibujar rectas por silicio ---
void DrawCalibLines(const CalData& pre, const CalData& post,
                    const char* title, const char* tagPre, const char* tagPost)
{
    int n = pre.a_vals.size();
    int nCols = 4;
    int nRows = (n + nCols - 1) / nCols;

    TCanvas* c = new TCanvas(title, title, 1600, 900);
    c->Divide(nCols, nRows);

    double chMin = 0, chMax = 4000;
    std::vector<double> x = {chMin, chMax};

    for (int i = 0; i < n; i++) {
        c->cd(i+1);

        std::vector<double> y_pre = {pre.a_vals[i]+pre.b_vals[i]*chMin,
                                     pre.a_vals[i]+pre.b_vals[i]*chMax};
        std::vector<double> y_post= {post.a_vals[i]+post.b_vals[i]*chMin,
                                     post.a_vals[i]+post.b_vals[i]*chMax};

        TGraph* gPre  = new TGraph(2, x.data(), y_pre.data());
        TGraph* gPost = new TGraph(2, x.data(), y_post.data());
        gPre->SetLineColor(kBlue); gPre->SetLineWidth(2);
        gPost->SetLineColor(kRed); gPost->SetLineWidth(2);

        gPre->SetTitle(TString::Format("%s - silicio %d;Canal;E [keV]", title, i));
        gPre->Draw("AL");
        gPost->Draw("L SAME");

        auto leg = new TLegend(0.55,0.75,0.9,0.9);
        leg->AddEntry(gPre, tagPre,"l");
        leg->AddEntry(gPost,tagPost,"l");
        leg->Draw();
    }
}

// --- Funcion principal ---
void compareCalibrations()
{
    // Leer capas
    CalData f0_pre = ReadCalibrationFileRaw("./Outputs/s2384_f0.dat");
    CalData f0_post = ReadCalibrationFileRaw("./Outputs/post_experiment_cals/s2384_f0_last.dat");

    CalData l0_pre = ReadCalibrationFileRaw("./Outputs/s2384_l0.dat");
    CalData l0_post = ReadCalibrationFileRaw("./Outputs/post_experiment_cals/s2384_l0_last.dat");

    CalData f2_pre = ReadCalibrationFileRaw("./Outputs/s2384_f2.dat");
    CalData f2_post = ReadCalibrationFileRaw("./Outputs/post_experiment_cals/s2384_f2_last.dat");

    CalData r0_pre = ReadCalibrationFileRaw("./Outputs/s2384_r0.dat");
    CalData r0_post = ReadCalibrationFileRaw("./Outputs/post_experiment_cals/s2384_r0_last.dat");

    // --- Comparación parámetros y canal 3000 ---
    DrawComparison(f0_pre, f0_post, "f0", "pre", "post");
    DrawComparison(l0_pre, l0_post, "l0", "pre", "post");
    DrawComparison(f2_pre, f2_post, "f2", "pre", "post");
    DrawComparison(r0_pre, r0_post, "r0", "pre", "post");

    // --- Comparación rectas ---
    DrawCalibLines(f0_pre, f0_post, "f0_lines", "pre", "post");
    DrawCalibLines(l0_pre, l0_post, "l0_lines", "pre", "post");
    DrawCalibLines(f2_pre, f2_post, "f2_lines", "pre", "post");
    DrawCalibLines(r0_pre, r0_post, "r0_lines", "pre", "post");
}
