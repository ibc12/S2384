#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TMultiGraph.h"

#include "FitModel.h"

#include <iostream>
#include <vector>

// Función auxiliar para evaluar Breit-Wigner con threshold
double BW(double x, double Er, double Gamma, double s)
{
    return Gamma * 0.159154943 / ((x - Er) * (x - Er) + Gamma * Gamma / 4);
}

double EvalManualConvolution(double x, double mean, double sigma, double Gamma0, double s)
{
    // Integration settings
    const int nPoints = 100;
    const double nSigmas = 4.0;

    double xMin = x - nSigmas * sigma;
    double xMax = x + nSigmas * sigma;
    double dx = (xMax - xMin) / nPoints;

    double result = 0.0;

    // Numerical integration of: BW(x') * Gauss(x - x') dx'
    for(int i = 0; i < nPoints; i++)
    {
        double xPrime = xMin + (i + 0.5) * dx;
        // Eval Breit-Wigner in xPrime (centered in mean)
        double bw = BW(xPrime, mean, Gamma0, s);
        // Eval Gaussian in (x - xPrime)
        double gauss = TMath::Gaus(x, xPrime, sigma, true); // normalized

        result += bw * gauss * dx;
    }
    return result;
}

void TestGammaConvolution()
{
    std::cout << "\n========================================" << std::endl;
    std::cout << "Testing Gamma(E) dependent Voigt" << std::endl;
    std::cout << "========================================\n" << std::endl;

    // Parámetros físicos para 7Li(d,p)
    double Sn = 2.03262;                                             // MeV - threshold energy
    double Er = 3.2;                                                 // MeV - resonance energy
    double Gamma0 = 1.0;                                             // MeV - resonance width at Er
    double sigma = 0.14;                                             // MeV - experimental resolution
    double mu = 7 * 1. / (7 + 1) * 931.5;                            // MeV/c^2 - reduced mass
    double R = (std::pow(7, 1. / 3.) + std::pow(1, 1. / 3.)) * 1.25; // fm - interaction radius
    int l = 2;                                                       // angular momentum

    std::cout << "Physical parameters:" << std::endl;
    std::cout << "  Sn (threshold) = " << Sn << " MeV" << std::endl;
    std::cout << "  Er (resonance) = " << Er << " MeV" << std::endl;
    std::cout << "  Gamma0         = " << Gamma0 << " MeV" << std::endl;
    std::cout << "  sigma          = " << sigma << " MeV" << std::endl;
    std::cout << "  l              = " << l << std::endl;
    std::cout << "  mu             = " << mu << " MeV/c^2" << std::endl;
    std::cout << "  R              = " << R << " fm" << std::endl;

    // Rango de energía para graficar
    double Emin = -4.0;
    double Emax = 8.0;
    int nPoints = 500;
    double dE = (Emax - Emin) / nPoints;

    // Vectores para almacenar datos
    std::vector<double> energy;
    std::vector<double> convolution; // Convolución final
    std::vector<double> voigt_root;  // TMath::Voigt para comparar

    std::cout << "\nCalculating functions..." << std::endl;

    for(int i = 0; i < nPoints; i++)
    {
        double E = Emin + (i + 0.5) * dE;
        energy.push_back(E);
        // 1. Convolution con Gamma constante
        double x = E;
        convolution.push_back(EvalManualConvolution(x, Er, sigma, Gamma0, l));
        // 2. TMath::Voigt para comparación (sin threshold)
        voigt_root.push_back(TMath::Voigt(E - Er, sigma, Gamma0));
    }

    std::cout << "Creating graphs..." << std::endl;

    // Crear TGraphs
    TGraph* gr_conv = new TGraph(nPoints, energy.data(), convolution.data());
    gr_conv->SetLineColor(kGreen + 2);
    gr_conv->SetLineWidth(3);
    gr_conv->SetTitle("Convolution (BW #otimes Gauss)");

    TGraph* gr_voigt = new TGraph(nPoints, energy.data(), voigt_root.data());
    gr_voigt->SetLineColor(kMagenta);
    gr_voigt->SetLineWidth(2);
    gr_voigt->SetLineStyle(3);
    gr_voigt->SetTitle("TMath::Voigt (no threshold)");

    // Canvas principal
    TCanvas* c1 = new TCanvas("c1", "Gamma(E) dependent Voigt Test", 1400, 800);
    c1->Divide(2, 2);

    // Panel 1: Todas las funciones
    c1->cd(1);
    gPad->SetGridy();
    TMultiGraph* mg1 = new TMultiGraph();
    mg1->Add(gr_conv, "L");
    mg1->Add(gr_voigt, "L");
    mg1->Draw("A");
    mg1->SetTitle("All functions;E_{x} (MeV);Amplitude (a.u.)");
    mg1->GetXaxis()->SetRangeUser(Emin, Emax);

    TLegend* leg1 = new TLegend(0.55, 0.55, 0.88, 0.88);
    leg1->AddEntry(gr_conv, "Convolution (final)", "l");
    leg1->AddEntry(gr_voigt, "ROOT Voigt (reference)", "l");
    leg1->Draw();

    // Línea vertical en el threshold
    TLine* line_sn = new TLine(Sn, 0, Sn, mg1->GetYaxis()->GetXmax() * 0.8);
    line_sn->SetLineColor(kBlack);
    line_sn->SetLineStyle(2);
    line_sn->SetLineWidth(2);
    line_sn->Draw();

    TLatex* tex_sn = new TLatex(Sn, mg1->GetYaxis()->GetXmax() * 0.85, "S_{n}");
    tex_sn->SetTextSize(0.04);
    tex_sn->Draw();

    // Panel 2: Zoom cerca del threshold
    c1->cd(2);
    gPad->SetGridy();
    TMultiGraph* mg2 = new TMultiGraph();
    mg2->Add((TGraph*)gr_conv->Clone(), "L");
    mg2->Draw("A");
    mg2->SetTitle("Zoom near threshold;E_{x} (MeV);Amplitude (a.u.)");
    mg2->GetXaxis()->SetRangeUser(Sn - 0.5, Sn + 1.0);

    TLegend* leg2 = new TLegend(0.55, 0.65, 0.88, 0.88);
    leg2->AddEntry(gr_conv, "Convolution (final)", "l");
    leg2->Draw();

    TLine* line_sn2 = new TLine(Sn, 0, Sn, mg2->GetYaxis()->GetXmax() * 0.9);
    line_sn2->SetLineColor(kBlack);
    line_sn2->SetLineStyle(2);
    line_sn2->SetLineWidth(2);
    line_sn2->Draw();

    // Panel 4: Ratio de convolución vs ROOT Voigt
    c1->cd(4);
    gPad->SetGridy();
    std::vector<double> ratio;
    for(size_t i = 0; i < convolution.size(); i++)
    {
        if(voigt_root[i] > 1e-10)
            ratio.push_back(convolution[i] / voigt_root[i]);
        else
            ratio.push_back(0);
    }

    TGraph* gr_ratio = new TGraph(nPoints, energy.data(), ratio.data());
    gr_ratio->SetLineColor(kBlack);
    gr_ratio->SetLineWidth(2);
    gr_ratio->Draw("AL");
    gr_ratio->SetTitle("Ratio: Convolution / ROOT Voigt;E_{x} (MeV);Ratio");
    gr_ratio->GetXaxis()->SetRangeUser(Emin, Emax);

    TLine* line_unity = new TLine(Emin, 1.0, Emax, 1.0);
    line_unity->SetLineColor(kRed);
    line_unity->SetLineStyle(2);
    line_unity->Draw();

    c1->Update();
    c1->SaveAs("test_gamma_convolution.pdf");

    std::cout << "\n========================================" << std::endl;
    std::cout << "Test completed successfully!" << std::endl;
    std::cout << "Output saved to: test_gamma_convolution.pdf" << std::endl;
    std::cout << "========================================\n" << std::endl;

}