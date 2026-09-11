#include "ActCutsManager.h"
#include "ActMergerData.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>

// FIX: Inclusión necesaria para controlar el directorio activo global de ROOT
#include <TH1D.h>
#include <TH2D.h>
#include <TPaveText.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "../../Fits/Histos.h"


// ============================================================================
// Obtiene el valor de straggling a partir del nombre del fichero
// ============================================================================
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

// ============================================================================
// Convierte la etiqueta del fichero en una cadena legible ("1-5" -> "1.5")
// ============================================================================
std::string FormatStragglingLabel(const std::string& label)
{
    std::string out = label;
    std::replace(out.begin(), out.end(), '-', '.');
    return out;
}

// ============================================================================
// Encuentra el índice del archivo de simulación según su etiqueta de straggling
// ============================================================================
int FindSimuFileIndex(const std::vector<std::string>& files, const std::string& targetLabel)
{
    for(size_t i = 0; i < files.size(); ++i)
    {
        auto label = FormatStragglingLabel(GetStragglingLabel(files[i]));
        if(label == targetLabel)
        {
            return static_cast<int>(i);
        }
    }
    return -1; // Retorna -1 si no se encuentra
}

// ============================================================================
// Divide un canvas en una rejilla aproximadamente cuadrada
// ============================================================================
void DivideCanvasGrid(TCanvas* c, int nPads)
{
    int nCols = static_cast<int>(std::ceil(std::sqrt(static_cast<double>(nPads))));
    int nRows = static_cast<int>(std::ceil(static_cast<double>(nPads) / nCols));
    c->Divide(nCols, nRows);
}

// ============================================================================
// Ajuste gaussiano
// ============================================================================
TF1* FitGaus1D(TH1D* h, double min, double max)
{
    static int fitCounter {0};
    auto name = "f_gaus_Ex_" + std::to_string(fitCounter++);

    auto* f = new TF1(name.c_str(), "gaus", min, max);

    // R = respeta el rango, Q = silencioso
    h->Fit(f, "RQ");

    f->SetLineColor(kRed);
    f->SetLineWidth(2);

    return f;
}

// ============================================================================
// Caja con los parámetros del ajuste
// ============================================================================
TPaveText* MakeParamBox(TF1* f)
{
    auto* pave = new TPaveText(0.60, 0.65, 0.89, 0.89, "NDC");

    pave->SetFillColor(0);
    pave->SetBorderSize(1);
    pave->SetTextFont(42);
    pave->SetTextSize(0.035);

    pave->AddText(TString::Format("Const = %.2f", f->GetParameter(0)));
    pave->AddText(TString::Format("Mean  = %.3f", f->GetParameter(1)));
    pave->AddText(TString::Format("Sigma = %.3f", f->GetParameter(2)));

    return pave;
}

// ============================================================================
// MAIN
// ============================================================================
void CompareResolutionsEx_L1_AngStr()
{
    // FIX: Evita que los histogramas creados posteriormente se asocien automáticamente a los archivos ROOT
    TH1::AddDirectory(kFALSE);

    // ========================================================================
    // PARTICULA
    // ========================================================================
    std::string particle = "p";

    // ========================================================================
    // FICHEROS DE SIMULACION
    // ========================================================================
    std::vector<std::string> files {};

    if(particle == "d")
    {
        files = {
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_1-5AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_2AngStr.root",
            //  "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_2-1AngStr.root",
            //  "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_2-2AngStr.root",
            //  "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_2-3AngStr.root",
            //  "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_2-4AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_2-5AngStr.root",
            // "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_2-6AngStr.root",
            // "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_2-7AngStr.root",
            // "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_2-8AngStr.root",
            // "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_2-9AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_3AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_3-1AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_3-2AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_3-3AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_3-4AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_3-5AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_4AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_4-5AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_5AngStr.root"};
    }
    else if(particle == "p")
    {
        files =
        {
            //"../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_1-5AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_2AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_2-5AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_2-7AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_3AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_3-3AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_3-5AngStr.root",
            
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_4-5AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_5AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_5-5AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_6AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_6-5AngStr.root",
            "../../Simulation/Outputs/7Li/test_ang_straggling_L1/2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_4AngStr.root"};
        }

        int wantedIdx = FindSimuFileIndex(files, "6.5");

        // ========================================================================
        // SIMULACION
        // ========================================================================
        std::vector<ROOT::RDataFrame> dfsSimu;
        dfsSimu.reserve(files.size());
        std::vector<ROOT::RDF::RResultPtr<TH1D>> histsExSimu;

        for(size_t i = 0; i < files.size(); i++)
        {
            auto label = FormatStragglingLabel(GetStragglingLabel(files[i]));
            dfsSimu.emplace_back(ROOT::RDataFrame {"SimulationTTree", files[i].c_str()});
            auto& df_simu = dfsSimu.back();
            auto nameEx = "hist_Ex_simu_" + std::to_string(i);
            auto titleEx = "Ex vs nPads (AngStr = " + label + "#circ); Ex [MeV]; nPads";
            histsExSimu.push_back(df_simu.Histo1D({nameEx.c_str(), titleEx.c_str(), 100, -2, 2}, "Eex"));
        }

        // ========================================================================
        // EXPERIMENTO
        // ========================================================================
        ROOT::RDataFrame df {
            "Final_Tree",
            TString::Format("../../PostAnalysis/Outputs/tree_ex_F_7Li_d_%s_filtered.root", particle.c_str())};
        auto def {df.Filter(
            [&particle](ActRoot::MergerData& m, double ex)
            {
                if(particle == "d")
                {
                    return m.fLight.IsFilled() == false;
                }
                else
                {
                    return m.fLight.IsFilled() == false;
                }
            },
            {"MergerData", "Ex"})};

        std::string nameEx = "hist_Ex_exp";
        std::string titleEx = "Ex vs nPads (Exp); Ex [MeV]; nPads";
        auto hEx = def.Histo1D({nameEx.c_str(), titleEx.c_str(), 50, -2, 2}, "Ex");

        // ========================================================================
        // EX 1D: SIMULACION + EXPERIMENTO
        // ========================================================================
        double exGausFitMin {-2.};
        double exGausFitMax {0.5};
        const int nTotalPads = static_cast<int>(histsExSimu.size()) + 1;
        std::vector<TF1*> fitsEx1D;
        std::vector<TPaveText*> pavesEx1D;
        std::vector<std::pair<std::string, double>> stragglingSigmas;

        // FIX: reset de directorio
        gROOT->cd();
        auto c_Ex1D = new TCanvas("c_Ex1D", "Ex (1D) - Ajuste gaussiano", 1600, 1000);
        DivideCanvasGrid(c_Ex1D, nTotalPads);

        for(size_t i = 0; i < histsExSimu.size(); i++)
        {
            c_Ex1D->cd(static_cast<int>(i) + 1);
            auto* h = static_cast<TH1D*>(histsExSimu[i]->DrawClone("HIST"));
            h->SetDirectory(0); // FIX: desvincular de ficheros
            auto* f = FitGaus1D(h, exGausFitMin, exGausFitMax);
            f->Draw("SAME");
            auto* pave = MakeParamBox(f);
            pave->Draw("SAME");
            fitsEx1D.push_back(f);
            pavesEx1D.push_back(pave);
            auto label = FormatStragglingLabel(GetStragglingLabel(files[i]));
            stragglingSigmas.emplace_back(label, f->GetParameter(2));
        }

        // ========================================================================
        // EXPERIMENTO COMPLETO
        // ========================================================================
        double sigmaExp {};
        c_Ex1D->cd(nTotalPads);
        {
            auto* h = static_cast<TH1D*>(hEx->DrawClone("HIST"));
            h->SetDirectory(0); // FIX
            auto* f = FitGaus1D(h, exGausFitMin, exGausFitMax);
            f->Draw("SAME");
            auto* pave = MakeParamBox(f);
            pave->Draw("SAME");
            fitsEx1D.push_back(f);
            pavesEx1D.push_back(pave);
            sigmaExp = f->GetParameter(2);
        }

        // ========================================================================
        // COMPARACION DE SIGMA SIMULACION VS EXPERIMENTO
        // ========================================================================
        std::sort(stragglingSigmas.begin(), stragglingSigmas.end(), [sigmaExp](const auto& a, const auto& b)
                  { return std::abs(a.second - sigmaExp) < std::abs(b.second - sigmaExp); });

        std::cout << "\n=== Comparacion sigma(Ex) simulacion vs experimento ===\n";
        std::cout << "Sigma experimental = " << sigmaExp << " MeV\n";

        for(const auto& [label, sigma] : stragglingSigmas)
        {
            std::cout << "  AngStr = " << label << "   sigma_simu = " << sigma
                      << " MeV   |dSigma| = " << std::abs(sigma - sigmaExp) << " MeV\n";
        }

        if(!stragglingSigmas.empty())
        {
            std::cout << "--> Mejor coincidencia: AngStr = " << stragglingSigmas.front().first << "\n";
        }

        // ========================================================================
        // EX EXPERIMENTAL POR FRANJAS DE TL
        // ========================================================================
        auto def_TL = def.Define("TL", [](ActRoot::MergerData& m) { return m.fLight.fTL; }, {"MergerData"});
        auto minTLRes = def_TL.Min<float>("TL");
        auto maxTLRes = def_TL.Max<float>("TL");
        double minTL = static_cast<double>(minTLRes.GetValue());
        double maxTL = static_cast<double>(maxTLRes.GetValue());
        const double TLbandWidth {15.};
        const int nTLbands = std::max(1, static_cast<int>(std::ceil((maxTL - minTL) / TLbandWidth)));
        std::vector<ROOT::RDF::RResultPtr<TH1D>> histsExByTL;
        std::vector<std::pair<double, double>> TLbandRanges;
        histsExByTL.reserve(nTLbands);
        TLbandRanges.reserve(nTLbands);

        for(int i = 0; i < nTLbands; i++)
        {
            double lowTL = minTL + i * TLbandWidth;
            double highTL = std::min(lowTL + TLbandWidth, maxTL);
            bool isLastBand = (i == nTLbands - 1);

            TLbandRanges.emplace_back(lowTL, highTL);

            auto defBand = def_TL.Filter(
                [lowTL, highTL, isLastBand](float TL)
                { return isLastBand ? (TL >= lowTL && TL <= highTL) : (TL >= lowTL && TL < highTL); }, {"TL"});

            auto nameExBand = "hist_Ex_TL_band_exp_" + std::to_string(i);
            auto titleExBand = TString::Format("Ex, TL #in [%.1f, %.1f) mm; Ex [MeV]; Counts", lowTL, highTL);

            histsExByTL.push_back(defBand.Histo1D({nameExBand.c_str(), titleExBand.Data(), 100, -2, 2}, "Ex"));
        }

        gROOT->cd(); // FIX: nos aseguramos de estar en la memoria raiz
        auto c_Ex_TL = new TCanvas("c_Ex_TL", "Ex experimental por franjas de TL - Ajustes gaussianos", 1600, 1000);
        DivideCanvasGrid(c_Ex_TL, nTLbands);

        std::vector<TF1*> fitsExTL;
        std::vector<TPaveText*> pavesExTL;

        std::cout << "\n=== Ex experimental por franjas de TL ===\n";
        std::cout << "Min TL = " << minTL << " mm, Max TL = " << maxTL << " mm, ancho = " << TLbandWidth << " mm\n";

        for(int i = 0; i < nTLbands; i++)
        {
            c_Ex_TL->cd(i + 1);

            auto* h = static_cast<TH1D*>(histsExByTL[i]->DrawClone("HIST"));
            h->SetDirectory(0); // FIX
            const auto& [lowTL, highTL] = TLbandRanges[i];

            if(h->GetEntries() < 10)
            {
                std::cout << "TL in [" << lowTL << ", " << highTL << ") mm : " << h->GetEntries()
                          << " entradas --> estadistica insuficiente, no se ajusta\n";
                fitsExTL.push_back(nullptr);
                continue;
            }

            auto* f = FitGaus1D(h, exGausFitMin, exGausFitMax);
            f->Draw("SAME");

            auto* pave = MakeParamBox(f);
            pave->Draw("SAME");
            fitsExTL.push_back(f);
            pavesExTL.push_back(pave);

            std::cout << "TL in [" << lowTL << ", " << highTL << ") mm : " << h->GetEntries() << " entradas"
                      << " --> Mean = " << f->GetParameter(1) << " +/- " << f->GetParError(1) << " MeV"
                      << ", Sigma = " << f->GetParameter(2) << " +/- " << f->GetParError(2) << " MeV\n";
        }

        // ========================================================================
        // TGraphErrors EXPERIMENTO
        // ========================================================================
        gROOT->cd(); // FIX
        auto c_Sigma_vs_TL = new TCanvas("c_Sigma_vs_TL", "Sigma(Ex) vs TL (Exp)", 800, 600);
        c_Sigma_vs_TL->cd();

        auto graph_Sigma_vs_TL = new TGraphErrors();
        graph_Sigma_vs_TL->SetName("graph_Sigma_vs_TL_Exp");
        graph_Sigma_vs_TL->SetTitle("Sigma(Ex) vs TL (Exp); TL [mm]; Sigma(Ex) [MeV]");
        graph_Sigma_vs_TL->SetMarkerStyle(20);
        graph_Sigma_vs_TL->SetMarkerSize(1.2);
        graph_Sigma_vs_TL->SetMarkerColor(kBlue);
        graph_Sigma_vs_TL->SetLineColor(kBlue);

        int pointIdxExp = 0;
        for(size_t i = 0; i < fitsExTL.size(); i++)
        {
            auto* f = fitsExTL[i];
            if(!f)
                continue;

            const auto& [lowTL, highTL] = TLbandRanges[i];
            double meanTL = (lowTL + highTL) / 2.0;
            double sigmaTL = (highTL - lowTL) / 2.0;

            double sigmaEx = f->GetParameter(2);
            double sigmaExErr = f->GetParError(2);

            graph_Sigma_vs_TL->SetPoint(pointIdxExp, meanTL, sigmaEx);
            graph_Sigma_vs_TL->SetPointError(pointIdxExp, sigmaTL, sigmaExErr);
            pointIdxExp++;
        }

        graph_Sigma_vs_TL->Draw("APE");

        // ========================================================================
        // ANALISIS TL SIMULACION (3-5 AngStr -> Índice 16)
        // ========================================================================
        // Seleccionamos el DataFrame correspondiente
        auto& df_simu_target = dfsSimu[wantedIdx];

        auto minTLSimu = df_simu_target.Min<double>("TL");
        auto maxTLSimu = df_simu_target.Max<double>("TL");
        double minTL_Simu = static_cast<double>(minTLSimu.GetValue());
        double maxTL_Simu = static_cast<double>(maxTLSimu.GetValue());
        const int nTLbands_Simu = std::max(1, static_cast<int>(std::ceil((maxTL_Simu - minTL_Simu) / TLbandWidth)));

        std::vector<ROOT::RDF::RResultPtr<TH1D>> histsExByTL_Simu;
        std::vector<std::pair<double, double>> TLbandRanges_Simu;
        histsExByTL_Simu.reserve(nTLbands_Simu);
        TLbandRanges_Simu.reserve(nTLbands_Simu);

        for(int i = 0; i < nTLbands_Simu; i++)
        {
            double lowTL = minTL_Simu + i * TLbandWidth;
            double highTL = std::min(lowTL + TLbandWidth, maxTL_Simu);
            bool isLastBand = (i == nTLbands_Simu - 1);

            TLbandRanges_Simu.emplace_back(lowTL, highTL);

            auto defBand_Simu = df_simu_target.Filter(
                [lowTL, highTL, isLastBand](double TL)
                { return isLastBand ? (TL >= lowTL && TL <= highTL) : (TL >= lowTL && TL < highTL); }, {"TL"});

            auto nameExSimuBand = "hist_Ex_TL_band_simu_" + std::to_string(i);
            auto titleExSimuBand =
                TString::Format("Ex, TL #in [%.1f, %.1f) mm (Simu 3-5 AngStr); Ex [MeV]; Counts", lowTL, highTL);

            histsExByTL_Simu.push_back(
                defBand_Simu.Histo1D({nameExSimuBand.c_str(), titleExSimuBand.Data(), 50, -2, 2}, "Eex"));
        }

        gROOT->cd(); // FIX
        auto c_Ex_TL_Simu = new TCanvas(
            "c_Ex_TL_Simu", "Ex simulation (3-5 AngStr) por franjas de TL - Ajustes gaussianos", 1600, 1000);
        DivideCanvasGrid(c_Ex_TL_Simu, nTLbands_Simu);

        std::vector<TF1*> fitsExTL_Simu;
        std::vector<TPaveText*> pavesExTL_Simu;

        for(int i = 0; i < nTLbands_Simu; i++)
        {
            c_Ex_TL_Simu->cd(i + 1);
            auto* h = static_cast<TH1D*>(histsExByTL_Simu[i]->DrawClone("HIST"));
            h->SetDirectory(0); // FIX: vital desvincular del TFile de la simulación
            const auto& [lowTL, highTL] = TLbandRanges_Simu[i];

            if(h->GetEntries() < 10)
            {
                std::cout << "Simu TL in [" << lowTL << ", " << highTL << ") mm : " << h->GetEntries()
                          << " entradas --> estadistica insuficiente, no se ajusta\n";
                fitsExTL_Simu.push_back(nullptr);
                continue;
            }

            auto* f = FitGaus1D(h, exGausFitMin, exGausFitMax);
            f->Draw("SAME");

            auto* pave = MakeParamBox(f);
            pave->Draw("SAME");
            fitsExTL_Simu.push_back(f);
            pavesExTL_Simu.push_back(pave);

            std::cout << "Simu TL in [" << lowTL << ", " << highTL << ") mm : " << h->GetEntries() << " entradas"
                      << " --> Mean = " << f->GetParameter(1) << " +/- " << f->GetParError(1) << " MeV"
                      << ", Sigma = " << f->GetParameter(2) << " +/- " << f->GetParError(2) << " MeV\n";
        }

        // ========================================================================
        // TGraphErrors SIMULACION
        // ========================================================================
        gROOT->cd(); // FIX: regresamos a la raíz global antes de crear el gráfico final
        auto c_Sigma_vs_TL_Simu = new TCanvas("c_Sigma_vs_TL_Simu", "Sigma(Ex) vs TL (Simu 3-5 AngStr)", 800, 600);
        c_Sigma_vs_TL_Simu->cd();

        auto graph_Sigma_vs_TL_Simu = new TGraphErrors();
        graph_Sigma_vs_TL_Simu->SetName("graph_Sigma_vs_TL_Simu");
        graph_Sigma_vs_TL_Simu->SetTitle("Sigma(Ex) vs TL (Simu 3-5 AngStr); TL [mm]; Sigma(Ex) [MeV]");
        graph_Sigma_vs_TL_Simu->SetMarkerStyle(21);
        graph_Sigma_vs_TL_Simu->SetMarkerSize(1.2);
        graph_Sigma_vs_TL_Simu->SetMarkerColor(kRed);
        graph_Sigma_vs_TL_Simu->SetLineColor(kRed);

        int pointIdxSimu = 0;
        for(size_t i = 0; i < fitsExTL_Simu.size(); i++)
        {
            auto* f = fitsExTL_Simu[i];
            if(!f)
                continue;

            const auto& [lowTL, highTL] = TLbandRanges_Simu[i];
            double meanTL = (lowTL + highTL) / 2.0;
            double sigmaTL = (highTL - lowTL) / 2.0;

            double sigmaEx = f->GetParameter(2);
            double sigmaExErr = f->GetParError(2);

            graph_Sigma_vs_TL_Simu->SetPoint(pointIdxSimu, meanTL, sigmaEx);
            graph_Sigma_vs_TL_Simu->SetPointError(pointIdxSimu, sigmaTL, sigmaExErr);
            pointIdxSimu++;
        }

        graph_Sigma_vs_TL_Simu->Draw("APE");
    }