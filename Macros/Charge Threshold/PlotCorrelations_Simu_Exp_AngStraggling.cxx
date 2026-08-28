#include "ActCutsManager.h"
#include "ActMergerData.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TPaveText.h>
#include <TProfile.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

constexpr int yMinExclusionZone = 56;
constexpr int yMaxExclusionZone = 71;

// Rangos de ajuste
constexpr double thetaFitMin = 76;
constexpr double thetaFitMax = 84;
constexpr double nPadsFitMin = 8;
constexpr double nPadsFitMax = 100;
constexpr double TLFitMin = 10;
constexpr double TLFitMax = 120;

// Rango de ajuste para la gaussiana del histograma 1D de theta.
// Por defecto igual al rango del ajuste pol1 (thetaFitMin/Max); ajustalo si
// el pico de theta no cae dentro de ese rango.
constexpr double thetaGausFitMin = 74;
constexpr double thetaGausFitMax = 86;

// Extrae una etiqueta corta (el token del angular straggling) de la ruta del
// fichero, p.ej. ".../L1_1-5AngStr.root" -> "1-5"
std::string GetStragglingLabel(const std::string& path)
{
    auto filename = path.substr(path.find_last_of('/') + 1);
    auto start = filename.find("_L1_");
    auto end = filename.find("AngStr");
    if(start != std::string::npos && end != std::string::npos)
    {
        start += 4; // saltar "_L1_"
        return filename.substr(start, end - start);
    }
    return filename;
}

// Convierte la etiqueta del fichero (que usa '-' en vez de '.' para el
// separador decimal, p.ej. "1-5") en una cadena legible para los titulos,
// p.ej. "1-5" -> "1.5"
std::string FormatStragglingLabel(const std::string& label)
{
    std::string out = label;
    std::replace(out.begin(), out.end(), '-', '.');
    return out;
}

// Divide un canvas en una rejilla lo mas cuadrada posible para nPads pads
void DivideCanvasGrid(TCanvas* c, int nPads)
{
    int cols = static_cast<int>(std::ceil(std::sqrt(nPads)));
    int rows = static_cast<int>(std::ceil(static_cast<double>(nPads) / cols));
    c->Divide(cols, rows);
}

// Hace el ProfileX (restringido en Y = nPads) y ajusta ese perfil con la
// funcion indicada ("pol1" o "pol2") en el rango [xMin, xMax].
// Devuelve el TF1 ya ajustado (con el resultado guardado dentro).
TF1* FitProfileX(TH2D* h, const std::string& fitFunc, double xMin, double xMax, double yMin, double yMax)
{
    int yBinLow = h->GetYaxis()->FindBin(yMin);
    int yBinHigh = h->GetYaxis()->FindBin(yMax);

    std::string profName = std::string(h->GetName()) + "_pfx";
    auto* prof = h->ProfileX(profName.c_str(), yBinLow, yBinHigh);
    prof->SetDirectory(nullptr);

    std::string fitName = std::string(h->GetName()) + "_fit";
    auto* f = new TF1(fitName.c_str(), fitFunc.c_str(), xMin, xMax);
    f->SetLineColor(kRed);
    f->SetLineWidth(2);

    // R: usa el rango del TF1; Q: silencioso; 0: no dibuja nada todavia
    prof->Fit(f, "RQ0");

    return f;
}

// Ajusta un histograma 1D con una gaussiana en el rango [xMin, xMax].
// Devuelve el TF1 ya ajustado (con el resultado guardado dentro).
TF1* FitGaus1D(TH1D* h, double xMin, double xMax)
{
    std::string fitName = std::string(h->GetName()) + "_fit";
    auto* f = new TF1(fitName.c_str(), "gaus", xMin, xMax);
    f->SetLineColor(kRed);
    f->SetLineWidth(2);

    // R: usa el rango del TF1; Q: silencioso; 0: no dibuja nada todavia
    h->Fit(f, "RQ0");

    return f;
}

// Construye una caja de texto con los parametros del ajuste y el chi2/ndf
TPaveText* MakeParamBox(TF1* f)
{
    auto* pave = new TPaveText(0.55, 0.65, 0.89, 0.89, "NDC");
    pave->SetFillColor(0);
    pave->SetBorderSize(1);
    pave->SetTextSize(0.045);
    pave->SetTextAlign(12);

    for(int i = 0; i < f->GetNpar(); i++)
        pave->AddText(Form("p%d = %.3g #pm %.3g", i, f->GetParameter(i), f->GetParError(i)));

    double ndf = f->GetNDF();
    pave->AddText(Form("#chi^{2}/ndf = %.2f", ndf > 0 ? f->GetChisquare() / ndf : 0.));

    return pave;
}

void PlotCorrelations_Simu_Exp_AngStraggling()
{
    std::string particle = "d";
    std::vector<std::string> files {};
    if(particle == "d")
    {
        // NOTA: revisa/ajusta esta lista para que coincida exactamente con
        // los ficheros que tengas en test_ang_straggling_L1 (nombres y
        // valores de straggling disponibles).
        files = {
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_1-5AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_2AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_2-1AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_2-2AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_2-3AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_2-4AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_2-5AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_2-6AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_2-7AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_2-8AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_2-9AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_3AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_3-1AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_3-2AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_3-3AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_3-4AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_3-5AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_4AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_4-5AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_5AngStr.root",
        };
    }
    else if(particle == "p")
    {
        // NOTA: revisa/ajusta esta lista para que coincida exactamente con
        // los ficheros que tengas en test_ang_straggling_L1 (nombres y
        // valores de straggling disponibles).
        files = {
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_1-5AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_2AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_2-5AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_3AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_3-5AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_4AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_4-5AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_5AngStr.root",
        };
    }

    // --- Datos experimentales (se construyen una sola vez) ---
    ROOT::DisableImplicitMT();
    ROOT::RDataFrame df_exp {
        "Final_Tree",
        TString::Format("../../PostAnalysis/Outputs/tree_ex_F_7Li_d_%s_filtered.root", particle.c_str()).Data()};
    auto df_exp_filtered =
        df_exp.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"});

    auto df_exp_nPads = df_exp_filtered
                            .Define("nPads",
                                    [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                                    {
                                        int nPads {};
                                        for(auto cluster : tpc.fClusters)
                                        {
                                            for(auto voxel : cluster.GetRefToVoxels())
                                            {
                                                int iy = voxel.GetPosition().Y();
                                                if(iy < yMinExclusionZone || iy > yMaxExclusionZone)
                                                    nPads++;
                                            }
                                        }
                                        return nPads;
                                    },
                                    {"MergerData", "TPCData"})
                            .Define("theta3Lab", [](ActRoot::MergerData& m) { return m.fThetaLight; }, {"MergerData"})
                            .Define("TL", [](ActRoot::MergerData& m) { return m.fLight.fTL; }, {"MergerData"});

    auto hist_theta_nPads_exp = df_exp_nPads.Histo2D(
        {"hist_theta_nPads_exp", "Theta vs nPads (Exp); theta3Lab [deg]; nPads", 400, 70, 190, 100, 0, 100}, "theta3Lab",
        "nPads");
    auto hist_TL_nPads_exp = df_exp_nPads.Histo2D(
        {"hist_TL_nPads_exp", "TL vs nPads (Exp); TL [mm]; nPads", 100, 0, 300, 100, 0, 100}, "TL", "nPads");
    auto hist_theta_exp =
        df_exp_nPads.Histo1D({"hist_theta_exp", "Theta (Exp); theta3Lab [deg]; Counts", 80, 70, 95}, "theta3Lab");

    // --- Bucle sobre todos los ficheros de simulacion ---
    std::vector<ROOT::RDataFrame> dfsSimu;
    dfsSimu.reserve(files.size());
    std::vector<ROOT::RDF::RResultPtr<TH2D>> histsThetaSimu;
    std::vector<ROOT::RDF::RResultPtr<TH2D>> histsTLSimu;
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histsTheta1DSimu;

    for(size_t i = 0; i < files.size(); i++)
    {
        auto label = FormatStragglingLabel(GetStragglingLabel(files[i]));

        dfsSimu.emplace_back(ROOT::RDataFrame {"SimulationTTree", files[i].c_str()});
        auto& df_simu = dfsSimu.back();

        auto nameTheta = "hist_theta_nPads_simu_" + std::to_string(i);
        auto titleTheta = "Theta vs nPads (AngStr = " + label + "#circ); theta3Lab [deg]; nPads";
        histsThetaSimu.push_back(
            df_simu.Histo2D({nameTheta.c_str(), titleTheta.c_str(), 400, 70, 190, 100, 0, 100}, "theta3Lab", "nPads"));

        auto nameTL = "hist_TL_nPads_simu_" + std::to_string(i);
        auto titleTL = "TL vs nPads (AngStr = " + label + "#circ); TL [mm]; nPads";
        histsTLSimu.push_back(
            df_simu.Histo2D({nameTL.c_str(), titleTL.c_str(), 100, 0, 300, 100, 0, 100}, "TL", "nPads"));

        auto nameTheta1D = "hist_theta1D_simu_" + std::to_string(i);
        auto titleTheta1D = "Theta (AngStr = " + label + "#circ); theta3Lab [deg]; Counts";
        histsTheta1DSimu.push_back(
            df_simu.Histo1D({nameTheta1D.c_str(), titleTheta1D.c_str(), 80, 70, 95}, "theta3Lab"));
    }

    int nTotalPads = static_cast<int>(files.size()) + 1; // +1 para el pad experimental

    // Vectores para mantener vivos los TF1 y TPaveText de cada pad
    std::vector<TF1*> fitsTheta, fitsTL;
    std::vector<TPaveText*> pavesTheta, pavesTL;

    // --- Canvas Theta vs nPads: un pad por straggling + uno para exp ---
    auto c_Theta_nPads = new TCanvas("c_Theta_nPads", "Theta vs nPads", 1600, 1000);
    DivideCanvasGrid(c_Theta_nPads, nTotalPads);
    for(size_t i = 0; i < histsThetaSimu.size(); i++)
    {
        c_Theta_nPads->cd(static_cast<int>(i) + 1);
        auto* h = histsThetaSimu[i]->DrawClone("COLZ");

        auto* f = FitProfileX(static_cast<TH2D*>(h), "pol1", thetaFitMin, thetaFitMax, nPadsFitMin, nPadsFitMax);
        f->Draw("SAME");
        auto* pave = MakeParamBox(f);
        pave->Draw("SAME");

        fitsTheta.push_back(f);
        pavesTheta.push_back(pave);
    }
    c_Theta_nPads->cd(nTotalPads);
    {
        auto* h = hist_theta_nPads_exp->DrawClone("COLZ");
        auto* f = FitProfileX(static_cast<TH2D*>(h), "pol1", thetaFitMin, thetaFitMax, nPadsFitMin, nPadsFitMax);
        f->Draw("SAME");
        auto* pave = MakeParamBox(f);
        pave->Draw("SAME");
        fitsTheta.push_back(f);
        pavesTheta.push_back(pave);
    }

    // --- Canvas TL vs nPads: un pad por straggling + uno para exp ---
    auto c_TL_nPads = new TCanvas("c_TL_nPads", "TL vs nPads", 1600, 1000);
    DivideCanvasGrid(c_TL_nPads, nTotalPads);
    for(size_t i = 0; i < histsTLSimu.size(); i++)
    {
        c_TL_nPads->cd(static_cast<int>(i) + 1);
        auto* h = histsTLSimu[i]->DrawClone("COLZ");

        auto* f = FitProfileX(static_cast<TH2D*>(h), "pol2", TLFitMin, TLFitMax, nPadsFitMin, nPadsFitMax);
        f->Draw("SAME");
        auto* pave = MakeParamBox(f);
        pave->Draw("SAME");

        fitsTL.push_back(f);
        pavesTL.push_back(pave);
    }
    c_TL_nPads->cd(nTotalPads);
    {
        auto* h = hist_TL_nPads_exp->DrawClone("COLZ");
        auto* f = FitProfileX(static_cast<TH2D*>(h), "pol2", TLFitMin, TLFitMax, nPadsFitMin, nPadsFitMax);
        f->Draw("SAME");
        auto* pave = MakeParamBox(f);
        pave->Draw("SAME");
        fitsTL.push_back(f);
        pavesTL.push_back(pave);
    }

    // --- Canvas Theta (1D): un pad por straggling + uno para exp, cada uno
    // ajustado con una gaussiana. Sirve para comparar visualmente y por
    // sigma que straggling se parece mas al experimento ---
    std::vector<TF1*> fitsTheta1D;
    std::vector<TPaveText*> pavesTheta1D;
    std::vector<std::pair<std::string, double>> stragglingSigmas; // (label, sigma)

    auto c_Theta1D = new TCanvas("c_Theta1D", "Theta (1D) - Ajuste gaussiano", 1600, 1000);
    DivideCanvasGrid(c_Theta1D, nTotalPads);
    for(size_t i = 0; i < histsTheta1DSimu.size(); i++)
    {
        c_Theta1D->cd(static_cast<int>(i) + 1);
        auto* h = histsTheta1DSimu[i]->DrawClone("HIST");

        auto* f = FitGaus1D(static_cast<TH1D*>(h), thetaGausFitMin, thetaGausFitMax);
        f->Draw("SAME");
        auto* pave = MakeParamBox(f);
        pave->Draw("SAME");

        fitsTheta1D.push_back(f);
        pavesTheta1D.push_back(pave);

        auto label = FormatStragglingLabel(GetStragglingLabel(files[i]));
        stragglingSigmas.emplace_back(label, f->GetParameter(2)); // p2 = sigma en "gaus"
    }

    double sigmaExp {};
    c_Theta1D->cd(nTotalPads);
    {
        auto* h = hist_theta_exp->DrawClone("HIST");
        auto* f = FitGaus1D(static_cast<TH1D*>(h), thetaGausFitMin, thetaGausFitMax);
        f->Draw("SAME");
        auto* pave = MakeParamBox(f);
        pave->Draw("SAME");
        fitsTheta1D.push_back(f);
        pavesTheta1D.push_back(pave);
        sigmaExp = f->GetParameter(2);
    }

    // Ordena las simulaciones por cercania de sigma a la experimental y
    // muestra un resumen para ver de un vistazo que straggling se parece mas
    std::sort(stragglingSigmas.begin(), stragglingSigmas.end(),
              [sigmaExp](const auto& a, const auto& b)
              { return std::abs(a.second - sigmaExp) < std::abs(b.second - sigmaExp); });

    std::cout << "\n=== Comparacion sigma(theta) simulacion vs experimento ===\n";
    std::cout << "Sigma experimental = " << sigmaExp << " deg\n";
    for(const auto& [label, sigma] : stragglingSigmas)
        std::cout << "  AngStr = " << label << "   sigma_simu = " << sigma
                   << " deg   |dSigma| = " << std::abs(sigma - sigmaExp) << " deg\n";
    if(!stragglingSigmas.empty())
        std::cout << "--> Mejor coincidencia: AngStr = " << stragglingSigmas.front().first << "\n";
}