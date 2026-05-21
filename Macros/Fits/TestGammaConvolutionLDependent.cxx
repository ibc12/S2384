#include "TCanvas.h"
#include "TF1.h"
#include "TF1Convolution.h" // NUEVO: Para usar TF1Convolution
#include "TGraph.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TMultiGraph.h"

#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

// ============================================================================
// Función para calcular Gamma(E) dependiente de l (penetrability factor)
// ============================================================================
std::function<double(double, double)> InitGammaL(int l, double s, double mu, double R)
{
    std::function<double(double, double)> ret;
    double hbar = 197.3269804; // MeV*fm

    if(l == 0)
    {
        ret = [s](double Ex, double Er)
        {
            double Ed = Ex - s;
            if(Ed <= 0)
                return 0.0; // Debajo del threshold
            return std::pow(Ed / Er, 0.5);
        };
    }
    else if(l == 1)
    {
        ret = [hbar, s, mu, R](double Ex, double Er)
        {
            double Ed = Ex - s;
            if(Ed <= 0)
                return 0.0; // Debajo del threshold

            return std::pow(Ed / Er, 1.5) * (2 * Ed / (Er + Ed)) * (1. + (2. * mu * Er * R * R) / (hbar * hbar)) /
                   (1. + (2. * mu * Ed * R * R) / (hbar * hbar));
        };
    }
    else if(l == 2)
    {
        ret = [hbar, s, mu, R](double Ex, double Er)
        {
            double Ed = Ex - s;
            if(Ed <= 0)
                return 0.0; // Debajo del threshold

            return std::pow(Ed / Er, 2.5) * (2 * Ed / (Er + Ed)) *
                   (9. + (6. * mu * Er * R * R) / (hbar * hbar) + std::pow((2. * mu * Er * R * R) / (hbar * hbar), 2)) /
                   (9. + (6. * mu * Ed * R * R) / (hbar * hbar) + std::pow((2. * mu * Ed * R * R) / (hbar * hbar), 2));
        };
    }
    else
    {
        throw std::runtime_error("InitGammaL(): only l = 0, 1, 2 implemented");
    }
    return ret;
}

// ============================================================================
// Breit-Wigner con threshold y Gamma(E)
// ============================================================================
double EvalBreitWigner(double Ex, double Er, double Gamma0, double s, std::function<double(double, double)>& gammaFunc)
{
    double Ed = Ex - s;

    // Si estamos debajo del threshold, BW = 0
    if(Ed <= 0)
        return 0.0;

    // Calcular Gamma dependiente de energía
    double Gamma = Gamma0 * gammaFunc(Ex, Er);

    // Breit-Wigner estándar (relativista)
    double denominator = std::pow(Ex - Er, 2) + std::pow(Gamma / 2.0, 2);
    return (Gamma / (2.0 * TMath::Pi())) / denominator;
}

// ============================================================================
// MÉTODO CORRECTO: Convolución manual BW(Gamma(E)) ⊗ Gauss
// ============================================================================
double EvalCorrectConvolution(double x, double Er, double sigma, double Gamma0, double s,
                              std::function<double(double, double)>& gammaFunc)
{
    // Integration settings
    const int nPoints = 15000;
    const double nSigmas = 45.0;

    double xMin = -2;
    double xMax = 9;
    double dx = 0.001;

    double result = 0.0;

    // Numerical integration: ∫ BW(x', Gamma(x')) * Gauss(x - x') dx'
    for(int i = 0; i < nPoints; i++)
    {
        double xPrime = xMin + (i + 0.5) * dx;

        // Evaluar Breit-Wigner con Gamma(xPrime) en xPrime
        double bw = EvalBreitWigner(xPrime, Er, Gamma0, s, gammaFunc);

        // Evaluar Gaussiana en (x - xPrime)
        double gauss = TMath::Gaus(x, xPrime, sigma, true); // normalized

        result += bw * gauss * dx;
    }
    return result;
}

// ============================================================================
// NUEVO: Convolución usando TF1Convolution de ROOT
// ============================================================================
double EvalRootConvolution(double x, double Er, double sigma, double Gamma0, double s,
                           std::function<double(double, double)>& gammaFunc)
{
    const double nSigmas = 4.0;

    // Define integration range
    double xMin = -2; // Asegurar que no vaya muy abajo del threshold
    double xMax = 9;

    // Crear lambda para Breit-Wigner compatible con TF1
    // Necesitamos capturar las variables por valor para que persistan
    auto bwLambda = [Er, Gamma0, s, &gammaFunc](double* xVal, double* /*par*/) -> double
    {
        double xPrime = xVal[0];
        return EvalBreitWigner(xPrime, Er, Gamma0, s, gammaFunc);
    };

    // Crear lambda para Gaussiana compatible con TF1
    // La gaussiana está centrada en 0 para la convolución
    auto gaussLambda = [sigma](double* xVal, double* /*par*/) -> double
    {
        double xPrime = xVal[0];
        return TMath::Gaus(xPrime, 0.0, sigma, true); // centered at 0, normalized
    };

    // Crear TF1 objects con nombres únicos para evitar conflictos
    static int call_counter = 0;
    call_counter++;

    TF1* fBW = new TF1(Form("fBW_temp_%d", call_counter), bwLambda, xMin, xMax, 0);
    TF1* fGauss = new TF1(Form("fGauss_temp_%d", call_counter), gaussLambda, xMin, xMax, 0);

    // Crear objeto de convolución
    TF1Convolution fConv(fBW, fGauss);

    // Configurar rango para la convolución
    fConv.SetRange(xMin, xMax);

    // Configurar número de puntos para FFT (más puntos = más precisión)
    fConv.SetNofPointsFFT(20000);

    // Crear TF1 desde la convolución
    TF1 fResult(Form("fResult_temp_%d", call_counter), fConv, xMin, xMax, 0);

    // Evaluar en el punto x
    double result = fResult.Eval(x);

    return result;
}

// ============================================================================
// MÉTODO INCORRECTO: Voigt con Gamma(x) modificado (lo que hacías antes)
// ============================================================================
double EvalIncorrectVoigt(double x, double Er, double sigma, double Gamma0, double s,
                          std::function<double(double, double)>& gammaFunc)
{
    // Calcular Gamma en el punto x (INCORRECTO)
    double gamma = Gamma0 * gammaFunc(x, Er);

    // Usar TMath::Voigt con ese Gamma (sin considerar threshold correctamente)
    return TMath::Voigt(x - Er, sigma, gamma);
}

// ============================================================================
// MAIN MACRO
// ============================================================================
void TestGammaConvolutionLDependent()
{
    std::cout << "\n========================================" << std::endl;
    std::cout << "Comparing CORRECT vs INCORRECT vs ROOT methods" << std::endl;
    std::cout << "for Gamma(E) dependent resonance (l=1)" << std::endl;
    std::cout << "========================================\n" << std::endl;

    // Parámetros físicos para 7Li(d,p)
    double Sn = 2.03262;                                             // MeV - threshold energy
    double Er = 2.2;                                                 // MeV - resonance energy (v1)
    double Gamma0 = 1;                                               // MeV - resonance width
    double sigma = 0.14;                                             // MeV - experimental resolution
    double mu = 7 * 1. / (7 + 1) * 931.5;                            // MeV/c^2 - reduced mass
    double R = (std::pow(7, 1. / 3.) + std::pow(1, 1. / 3.)) * 1.25; // fm - interaction radius
    int l = 1;                                                       // angular momentum

    std::cout << "Physical parameters:" << std::endl;
    std::cout << "  Sn (threshold) = " << Sn << " MeV" << std::endl;
    std::cout << "  Er (resonance) = " << Er << " MeV" << std::endl;
    std::cout << "  Gamma0         = " << Gamma0 << " MeV" << std::endl;
    std::cout << "  sigma          = " << sigma << " MeV" << std::endl;
    std::cout << "  l              = " << l << std::endl;
    std::cout << "  mu             = " << mu << " MeV/c^2" << std::endl;
    std::cout << "  R              = " << R << " fm" << std::endl;

    // Crear función Gamma(E) para l=1
    auto gammaFunc = InitGammaL(l, Sn, mu, R);

    // Rango de energía para graficar
    double Emin = 1.0;
    double Emax = 7.0;
    int nPoints = 100;
    double dE = (Emax - Emin) / nPoints;

    // Vectores para almacenar datos
    std::vector<double> energy;
    std::vector<double> gamma_values;     // Valores de Gamma(E)
    std::vector<double> bw_gamma_e;       // BW con Gamma(E) sin convolución
    std::vector<double> correct_method;   // CORRECTO: BW(Gamma(E)) ⊗ Gauss (manual)
    std::vector<double> root_method;      // NUEVO: BW(Gamma(E)) ⊗ Gauss (ROOT TF1Convolution)
    std::vector<double> incorrect_method; // INCORRECTO: Voigt(Gamma(E))
    std::vector<double> voigt_constant;   // Referencia: Voigt con Gamma constante

    std::cout << "\nCalculating functions..." << std::endl;

    for(int i = 0; i < nPoints; i++)
    {
        double E = Emin + (i + 0.5) * dE;
        energy.push_back(E);

        // 1. Gamma(E)
        double GammaE = Gamma0 * gammaFunc(E, Er);
        gamma_values.push_back(GammaE);

        // 2. BW con Gamma(E) sin convolución
        bw_gamma_e.push_back(EvalBreitWigner(E, Er, Gamma0, Sn, gammaFunc));

        // 3. MÉTODO CORRECTO: Convolución manual BW(Gamma(E)) ⊗ Gauss
        correct_method.push_back(EvalCorrectConvolution(E, Er, sigma, Gamma0, Sn, gammaFunc));

        // 4. NUEVO: Convolución usando ROOT TF1Convolution
        root_method.push_back(EvalRootConvolution(E, Er, sigma, Gamma0, Sn, gammaFunc));

        // 5. MÉTODO INCORRECTO: Voigt con Gamma(E) modificado
        incorrect_method.push_back(EvalIncorrectVoigt(E, Er, sigma, Gamma0, Sn, gammaFunc));

        // 6. Voigt con Gamma constante (referencia)
        voigt_constant.push_back(TMath::Voigt(E - Er, sigma, Gamma0));

        if((i + 1) % 20 == 0)
            std::cout << "  Progress: " << (i + 1) << "/" << nPoints << std::endl;
    }

    std::cout << "Creating graphs..." << std::endl;

    // Crear TGraphs
    TGraph* gr_gamma = new TGraph(nPoints, energy.data(), gamma_values.data());
    gr_gamma->SetLineColor(kOrange + 1);
    gr_gamma->SetLineWidth(2);
    gr_gamma->SetTitle("#Gamma(E)");

    TGraph* gr_bw = new TGraph(nPoints, energy.data(), bw_gamma_e.data());
    gr_bw->SetLineColor(kBlue);
    gr_bw->SetLineWidth(2);
    gr_bw->SetLineStyle(2);
    gr_bw->SetTitle("BW(#Gamma(E)) - no convolution");

    TGraph* gr_correct = new TGraph(nPoints, energy.data(), correct_method.data());
    gr_correct->SetLineColor(kGreen + 2);
    gr_correct->SetLineWidth(3);
    gr_correct->SetTitle("CORRECT: BW(#Gamma(E)) #otimes Gauss (manual)");

    TGraph* gr_root = new TGraph(nPoints, energy.data(), root_method.data());
    gr_root->SetLineColor(kCyan + 2);
    gr_root->SetLineWidth(2);
    gr_root->SetLineStyle(1);
    gr_root->SetTitle("ROOT: BW(#Gamma(E)) #otimes Gauss (TF1Convolution)");

    TGraph* gr_incorrect = new TGraph(nPoints, energy.data(), incorrect_method.data());
    gr_incorrect->SetLineColor(kRed);
    gr_incorrect->SetLineWidth(2);
    gr_incorrect->SetLineStyle(7);
    gr_incorrect->SetTitle("INCORRECT: Voigt(#Gamma(E))");

    TGraph* gr_voigt_const = new TGraph(nPoints, energy.data(), voigt_constant.data());
    gr_voigt_const->SetLineColor(kMagenta);
    gr_voigt_const->SetLineWidth(2);
    gr_voigt_const->SetLineStyle(3);
    gr_voigt_const->SetTitle("Reference: Voigt(#Gamma_{0})");

    // Canvas principal con más paneles
    TCanvas* c1 = new TCanvas("c1", "Gamma(E) Methods Comparison", 1800, 1200);
    c1->Divide(3, 3);

    // ========================================================================
    // Panel 1: Gamma(E) vs E
    // ========================================================================
    c1->cd(1);
    gPad->SetGridy();
    gr_gamma->Draw("AL");
    gr_gamma->GetXaxis()->SetTitle("E_{x} (MeV)");
    gr_gamma->GetYaxis()->SetTitle("#Gamma(E) (MeV)");
    gr_gamma->SetTitle("#Gamma(E) for l=1");
    gr_gamma->GetXaxis()->SetRangeUser(Emin, Emax);

    TLine* line_sn1 = new TLine(Sn, 0, Sn, *std::max_element(gamma_values.begin(), gamma_values.end()) * 1.1);
    line_sn1->SetLineColor(kBlack);
    line_sn1->SetLineStyle(2);
    line_sn1->SetLineWidth(2);
    line_sn1->Draw();

    TLatex* tex1 = new TLatex(0.15, 0.85, Form("l = %d, S_{n} = %.3f MeV", l, Sn));
    tex1->SetNDC();
    tex1->SetTextSize(0.04);
    tex1->Draw();

    // ========================================================================
    // Panel 2: BW con Gamma(E) (sin convolución)
    // ========================================================================
    c1->cd(2);
    gPad->SetGridy();
    gr_bw->Draw("AL");
    gr_bw->GetXaxis()->SetTitle("E_{x} (MeV)");
    gr_bw->GetYaxis()->SetTitle("Amplitude (a.u.)");
    gr_bw->SetTitle("BW with #Gamma(E) - NO convolution");
    gr_bw->GetXaxis()->SetRangeUser(Emin, Emax);

    TLine* line_sn2 = new TLine(Sn, 0, Sn, gr_bw->GetYaxis()->GetXmax() * 0.9);
    line_sn2->SetLineColor(kBlack);
    line_sn2->SetLineStyle(2);
    line_sn2->SetLineWidth(2);
    line_sn2->Draw();

    // ========================================================================
    // Panel 3: Comparación TODOS los métodos
    // ========================================================================
    c1->cd(3);
    gPad->SetGridy();
    TMultiGraph* mg_all = new TMultiGraph();
    mg_all->Add(gr_correct, "L");
    mg_all->Add(gr_root, "L");
    mg_all->Add(gr_incorrect, "L");
    mg_all->Add(gr_voigt_const, "L");
    mg_all->Draw("A");
    mg_all->SetTitle("All Methods Comparison;E_{x} (MeV);Amplitude (a.u.)");
    mg_all->GetXaxis()->SetRangeUser(Emin, Emax);

    TLegend* leg3 = new TLegend(0.40, 0.55, 0.88, 0.88);
    leg3->SetTextSize(0.030);
    leg3->AddEntry(gr_correct, "Manual convolution", "l");
    leg3->AddEntry(gr_root, "ROOT TF1Convolution", "l");
    leg3->AddEntry(gr_incorrect, "INCORRECT method", "l");
    leg3->AddEntry(gr_voigt_const, "Voigt(#Gamma_{0})", "l");
    leg3->Draw();

    TLine* line_sn3 = new TLine(Sn, 0, Sn, mg_all->GetYaxis()->GetXmax() * 0.9);
    line_sn3->SetLineColor(kBlack);
    line_sn3->SetLineStyle(2);
    line_sn3->SetLineWidth(2);
    line_sn3->Draw();

    // ========================================================================
    // Panel 4: Zoom cerca del threshold - Comparación Manual vs ROOT
    // ========================================================================
    c1->cd(4);
    gPad->SetGridy();
    TMultiGraph* mg_zoom = new TMultiGraph();
    mg_zoom->Add((TGraph*)gr_correct->Clone(), "L");
    mg_zoom->Add((TGraph*)gr_root->Clone(), "L");
    mg_zoom->Add((TGraph*)gr_bw->Clone(), "L");
    mg_zoom->Draw("A");
    mg_zoom->SetTitle("Zoom near threshold;E_{x} (MeV);Amplitude (a.u.)");
    mg_zoom->GetXaxis()->SetRangeUser(Sn - 0.3, Sn + 0.8);

    TLegend* leg4 = new TLegend(0.45, 0.65, 0.88, 0.88);
    leg4->SetTextSize(0.035);
    leg4->AddEntry(gr_correct, "Manual", "l");
    leg4->AddEntry(gr_root, "ROOT", "l");
    leg4->AddEntry(gr_bw, "BW (no conv)", "l");
    leg4->Draw();

    TLine* line_sn4 = new TLine(Sn, 0, Sn, mg_zoom->GetYaxis()->GetXmax());
    line_sn4->SetLineColor(kBlack);
    line_sn4->SetLineStyle(2);
    line_sn4->SetLineWidth(2);
    line_sn4->Draw();

    // ========================================================================
    // Panel 5: Ratio Manual / ROOT (verificación)
    // ========================================================================
    c1->cd(5);
    gPad->SetGridy();
    std::vector<double> ratio_manual_root;
    for(size_t i = 0; i < correct_method.size(); i++)
    {
        if(root_method[i] > 1e-10)
            ratio_manual_root.push_back(correct_method[i] / root_method[i]);
        else
            ratio_manual_root.push_back(1.0);
    }

    TGraph* gr_ratio_mr = new TGraph(nPoints, energy.data(), ratio_manual_root.data());
    gr_ratio_mr->SetLineColor(kBlue);
    gr_ratio_mr->SetLineWidth(2);
    gr_ratio_mr->Draw("AL");
    gr_ratio_mr->SetTitle("Ratio: Manual / ROOT;E_{x} (MeV);Ratio");
    gr_ratio_mr->GetXaxis()->SetRangeUser(Emin, Emax);

    TLine* line_unity1 = new TLine(Emin, 1.0, Emax, 1.0);
    line_unity1->SetLineColor(kRed);
    line_unity1->SetLineStyle(2);
    line_unity1->Draw();

    TLine* line_sn5 = new TLine(Sn, gr_ratio_mr->GetYaxis()->GetXmin(), Sn, gr_ratio_mr->GetYaxis()->GetXmax());
    line_sn5->SetLineColor(kBlack);
    line_sn5->SetLineStyle(2);
    line_sn5->SetLineWidth(2);
    line_sn5->Draw();

    // ========================================================================
    // Panel 6: Ratio CORRECT / INCORRECT
    // ========================================================================
    c1->cd(6);
    gPad->SetGridy();
    std::vector<double> ratio_methods;
    for(size_t i = 0; i < correct_method.size(); i++)
    {
        if(incorrect_method[i] > 1e-10)
            ratio_methods.push_back(correct_method[i] / incorrect_method[i]);
        else
            ratio_methods.push_back(1.0);
    }

    TGraph* gr_ratio = new TGraph(nPoints, energy.data(), ratio_methods.data());
    gr_ratio->SetLineColor(kBlack);
    gr_ratio->SetLineWidth(2);
    gr_ratio->Draw("AL");
    gr_ratio->SetTitle("Ratio: CORRECT / INCORRECT;E_{x} (MeV);Ratio");
    gr_ratio->GetXaxis()->SetRangeUser(Emin, Emax);

    TLine* line_unity2 = new TLine(Emin, 1.0, Emax, 1.0);
    line_unity2->SetLineColor(kRed);
    line_unity2->SetLineStyle(2);
    line_unity2->Draw();

    TLine* line_sn6 = new TLine(Sn, gr_ratio->GetYaxis()->GetXmin(), Sn, gr_ratio->GetYaxis()->GetXmax());
    line_sn6->SetLineColor(kBlack);
    line_sn6->SetLineStyle(2);
    line_sn6->SetLineWidth(2);
    line_sn6->Draw();

    // ========================================================================
    // Panel 7: Diferencia absoluta Manual - ROOT
    // ========================================================================
    c1->cd(7);
    gPad->SetGridy();
    std::vector<double> diff_manual_root;
    for(size_t i = 0; i < correct_method.size(); i++)
    {
        diff_manual_root.push_back(correct_method[i] - root_method[i]);
    }

    TGraph* gr_diff_mr = new TGraph(nPoints, energy.data(), diff_manual_root.data());
    gr_diff_mr->SetLineColor(kCyan + 2);
    gr_diff_mr->SetLineWidth(2);
    gr_diff_mr->Draw("AL");
    gr_diff_mr->SetTitle("Difference: Manual - ROOT;E_{x} (MeV);#Delta Amplitude");
    gr_diff_mr->GetXaxis()->SetRangeUser(Emin, Emax);

    TLine* line_zero1 = new TLine(Emin, 0, Emax, 0);
    line_zero1->SetLineColor(kRed);
    line_zero1->SetLineStyle(2);
    line_zero1->Draw();

    TLine* line_sn7 = new TLine(Sn, gr_diff_mr->GetYaxis()->GetXmin(), Sn, gr_diff_mr->GetYaxis()->GetXmax());
    line_sn7->SetLineColor(kBlack);
    line_sn7->SetLineStyle(2);
    line_sn7->SetLineWidth(2);
    line_sn7->Draw();

    // ========================================================================
    // Panel 8: Diferencia absoluta CORRECT - INCORRECT
    // ========================================================================
    c1->cd(8);
    gPad->SetGridy();
    std::vector<double> diff_absolute;
    for(size_t i = 0; i < correct_method.size(); i++)
    {
        diff_absolute.push_back(correct_method[i] - incorrect_method[i]);
    }

    TGraph* gr_diff = new TGraph(nPoints, energy.data(), diff_absolute.data());
    gr_diff->SetLineColor(kBlue);
    gr_diff->SetLineWidth(2);
    gr_diff->Draw("AL");
    gr_diff->SetTitle("Difference: CORRECT - INCORRECT;E_{x} (MeV);#Delta Amplitude");
    gr_diff->GetXaxis()->SetRangeUser(Emin, Emax);

    TLine* line_zero2 = new TLine(Emin, 0, Emax, 0);
    line_zero2->SetLineColor(kRed);
    line_zero2->SetLineStyle(2);
    line_zero2->Draw();

    TLine* line_sn8 = new TLine(Sn, gr_diff->GetYaxis()->GetXmin(), Sn, gr_diff->GetYaxis()->GetXmax());
    line_sn8->SetLineColor(kBlack);
    line_sn8->SetLineStyle(2);
    line_sn8->SetLineWidth(2);
    line_sn8->Draw();

    // ========================================================================
    // Panel 9: Zoom en el pico - Manual vs ROOT vs INCORRECT
    // ========================================================================
    c1->cd(9);
    gPad->SetGridy();
    TMultiGraph* mg_peak = new TMultiGraph();
    mg_peak->Add((TGraph*)gr_correct->Clone(), "L");
    mg_peak->Add((TGraph*)gr_root->Clone(), "L");
    mg_peak->Add((TGraph*)gr_incorrect->Clone(), "L");
    mg_peak->Draw("A");
    mg_peak->SetTitle("Zoom at resonance peak;E_{x} (MeV);Amplitude (a.u.)");
    mg_peak->GetXaxis()->SetRangeUser(Er - 0.5, Er + 0.5);

    TLegend* leg9 = new TLegend(0.45, 0.70, 0.88, 0.88);
    leg9->SetTextSize(0.035);
    leg9->AddEntry(gr_correct, "Manual", "l");
    leg9->AddEntry(gr_root, "ROOT", "l");
    leg9->AddEntry(gr_incorrect, "INCORRECT", "l");
    leg9->Draw();

    c1->Update();
    c1->SaveAs("comparison_all_methods.pdf");

    std::cout << "\n========================================" << std::endl;
    std::cout << "Comparison completed!" << std::endl;
    std::cout << "Output: comparison_all_methods.pdf" << std::endl;
    std::cout << "========================================\n" << std::endl;

    // ========================================================================
    // Print numerical comparison
    // ========================================================================
    std::cout << "\nNumerical comparison at key energies:" << std::endl;
    std::cout << "E (MeV) | Manual    | ROOT      | INCORRECT | M/R Ratio | Diff(M-R)% | M/I Ratio | Diff(M-I)%"
              << std::endl;
    std::cout << "--------|-----------|-----------|-----------|-----------|------------|-----------|------------"
              << std::endl;

    std::vector<double> test_energies = {Sn - 0.2, Sn, Sn + 0.3, Er - 0.5, Er, Er + 0.5, Er + 1.0};
    for(double E : test_energies)
    {
        int idx = (E - Emin) / dE;
        if(idx >= 0 && idx < nPoints)
        {
            double ratio_mr = (root_method[idx] > 1e-10) ? correct_method[idx] / root_method[idx] : 1.0;
            double diff_mr_pct =
                (root_method[idx] > 1e-10) ? 100.0 * (correct_method[idx] - root_method[idx]) / root_method[idx] : 0.0;
            double ratio_mi = (incorrect_method[idx] > 1e-10) ? correct_method[idx] / incorrect_method[idx] : 1.0;
            double diff_mi_pct = (incorrect_method[idx] > 1e-10)
                                     ? 100.0 * (correct_method[idx] - incorrect_method[idx]) / incorrect_method[idx]
                                     : 0.0;
            printf("%6.3f  | %9.6f | %9.6f | %9.6f | %9.6f | %9.4f%% | %9.6f | %9.4f%%\n", E, correct_method[idx],
                   root_method[idx], incorrect_method[idx], ratio_mr, diff_mr_pct, ratio_mi, diff_mi_pct);
        }
    }
}