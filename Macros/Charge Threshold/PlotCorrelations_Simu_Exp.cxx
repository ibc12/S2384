#include "ActMergerData.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"

#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TF1.h>
#include <TPaveText.h>

#include <cmath>
#include <iostream>
#include <string>
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

// Extrae una etiqueta corta (el token del threshold) de la ruta del fichero,
// p.ej. ".../L1_9e5Thresh.root" -> "9e5"
std::string GetThresholdLabel(const std::string& path)
{
    auto filename = path.substr(path.find_last_of('/') + 1);
    auto start = filename.find("_L1_");
    auto end = filename.find("Thresh");
    if(start != std::string::npos && end != std::string::npos)
    {
        start += 4; // saltar "_L1_"
        return filename.substr(start, end - start);
    }
    return filename;
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

void PlotCorrelations_Simu_Exp()
{
    std::vector<std::string> files = {
        "../../Simulation/Outputs/7Li/test_charge_threshold/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_9e5Thresh.root",
        "../../Simulation/Outputs/7Li/test_charge_threshold/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_8e5Thresh.root",
        "../../Simulation/Outputs/7Li/test_charge_threshold/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_7e5Thresh.root",
        "../../Simulation/Outputs/7Li/test_charge_threshold/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_6e5Thresh.root",
        "../../Simulation/Outputs/7Li/test_charge_threshold/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_5e5Thresh.root",
        "../../Simulation/Outputs/7Li/test_charge_threshold/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_4e5Thresh.root",
        "../../Simulation/Outputs/7Li/test_charge_threshold/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_3e5Thresh.root",
        "../../Simulation/Outputs/7Li/test_charge_threshold/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_2e5Thresh.root",
        "../../Simulation/Outputs/7Li/test_charge_threshold/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_1e5Thresh.root",
        "../../Simulation/Outputs/7Li/test_charge_threshold/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_9e4Thresh.root",
        "../../Simulation/Outputs/7Li/test_charge_threshold/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_8e4Thresh.root",
        "../../Simulation/Outputs/7Li/test_charge_threshold/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_7e4Thresh.root",
        "../../Simulation/Outputs/7Li/test_charge_threshold/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_6e4Thresh.root",
        "../../Simulation/Outputs/7Li/test_charge_threshold/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_5e4Thresh.root",
        "../../Simulation/Outputs/7Li/test_charge_threshold/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_4e4Thresh.root",
        "../../Simulation/Outputs/7Li/test_charge_threshold/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_3e4Thresh.root",
        "../../Simulation/Outputs/7Li/test_charge_threshold/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_2e4Thresh.root",
        "../../Simulation/Outputs/7Li/test_charge_threshold/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_1e4Thresh.root",
    };

    // --- Datos experimentales (se construyen una sola vez) ---
    ROOT::RDataFrame df_exp {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_7Li_d_d_filtered.root"};
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
        {"hist_theta_nPads_exp", "Theta vs nPads (Exp); theta3Lab [deg]; nPads", 100, 70, 90, 100, 0, 100}, "theta3Lab",
        "nPads");
    auto hist_TL_nPads_exp = df_exp_nPads.Histo2D(
        {"hist_TL_nPads_exp", "TL vs nPads (Exp); TL [mm]; nPads", 100, 0, 300, 100, 0, 100}, "TL", "nPads");

    // --- Bucle sobre todos los ficheros de simulacion ---
    std::vector<ROOT::RDataFrame> dfsSimu;
    dfsSimu.reserve(files.size());
    std::vector<ROOT::RDF::RResultPtr<TH2D>> histsThetaSimu;
    std::vector<ROOT::RDF::RResultPtr<TH2D>> histsTLSimu;

    for(size_t i = 0; i < files.size(); i++)
    {
        auto label = GetThresholdLabel(files[i]);

        dfsSimu.emplace_back(ROOT::RDataFrame {"SimulationTTree", files[i].c_str()});
        auto& df_simu = dfsSimu.back();

        auto nameTheta = "hist_theta_nPads_simu_" + std::to_string(i);
        auto titleTheta = "Theta vs nPads (Simu " + label + "); theta3Lab [deg]; nPads";
        histsThetaSimu.push_back(
            df_simu.Histo2D({nameTheta.c_str(), titleTheta.c_str(), 100, 70, 90, 100, 0, 100}, "theta3Lab", "nPads"));

        auto nameTL = "hist_TL_nPads_simu_" + std::to_string(i);
        auto titleTL = "TL vs nPads (Simu " + label + "); TL [mm]; nPads";
        histsTLSimu.push_back(
            df_simu.Histo2D({nameTL.c_str(), titleTL.c_str(), 100, 0, 300, 100, 0, 100}, "TL", "nPads"));
    }

    int nTotalPads = static_cast<int>(files.size()) + 1; // +1 para el pad experimental

    // Vectores para mantener vivos los TF1 y TPaveText de cada pad
    std::vector<TF1*> fitsTheta, fitsTL;
    std::vector<TPaveText*> pavesTheta, pavesTL;

    // --- Canvas Theta vs nPads: un pad por threshold + uno para exp ---
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

    // --- Canvas TL vs nPads: un pad por threshold + uno para exp ---
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
}