// ---------------------------------------------------------------------
// Fit of a relative-energy (Enf) TH1 histogram with:
//   (a) the s-wave function (Eqs. 1-2 of Johansson et al., NPA 842 (2010)
//       15-32) convoluted with the experimental resolution (Gaussian,
//       fixed sigma), and
//   (b) a Phase Space (PS) shape taken from a simulated histogram,
//       added with its own free amplitude.
//
// Fit parameters:
//   [0] Amplitude_sig -> FREE  (normalization of the s-wave component)
//   [1] as (fm)       -> FREE  (the scattering length to be extracted)
//   [2] r0 (fm)       -> FIXED
//   [3] eps (MeV)     -> FIXED
//   [4] sigma_res     -> FIXED (experimental resolution, MeV)
//   [5] Amplitude_PS  -> FREE  (normalization of the phase-space component)
//
// Usage in ROOT:
//   root -l fit_spectra_12Li_sWave.C
// (or wrap it as a function call with arguments if you prefer)
// ---------------------------------------------------------------------

#include "ROOT/RDataFrame.hxx"

#include <TCanvas.h>
#include <TF1.h>
#include <TF1Convolution.h>
#include <TFile.h>
#include <TFitResultPtr.h>
#include <TH1.h>
#include <TLegend.h>
#include <TMath.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <iostream>

#include "../../Fits/Histos.h"

static const double HBAR_C = 197.3269804; // MeV*fm

double MuC2_11Li_n()
{
    const double m_n_u = 1.008665;
    const double m_11Li_u = 11.042078;
    const double u_to_MeV = 931.494;
    return (m_n_u * m_11Li_u) / (m_n_u + m_11Li_u) * u_to_MeV;
}
static const double MU_C2 = MuC2_11Li_n();

double WaveNumber(double E)
{
    return TMath::Sqrt(2.0 * MU_C2 * E) / HBAR_C;
}

// Physical function (Eqs. 1-2), returns 0 for Enf <= 0.
// Parameters: [0]=as, [1]=r0, [2]=eps
double SWaveSpectrum(double* x, double* par)
{
    double Enf = x[0];
    if(Enf <= 0.0)
        return 0.0;

    double a_s = par[0];
    double r0 = par[1];
    double eps = par[2];

    double knf = WaveNumber(Enf);
    if(knf < 1e-8)
        knf = 1e-8;

    double k = WaveNumber(eps);

    double cot_delta = (-1.0 / a_s + 0.5 * r0 * knf * knf) / knf;
    double sin_delta = 1.0 / TMath::Sqrt(1.0 + cot_delta * cot_delta);
    double cos_delta = cot_delta / TMath::Sqrt(1.0 + cot_delta * cot_delta);

    double bracket = (1.0 / (k * k + knf * knf)) * cos_delta + (k / knf) * sin_delta;

    return knf * bracket * bracket;
}

// --- Global objects (must remain alive during the fit) ---
TF1* g_f_orig = nullptr;
TF1* g_f_gaus = nullptr;
TF1Convolution* g_conv = nullptr;
TF1* g_f_conv = nullptr; // s-wave convolution, without amplitude normalization
TH1* g_hPS = nullptr;    // phase-space shape (normalized to unit integral), no TF1 wrapper needed

// Evaluate the PS shape at x, returning 0 outside the histogram range.
// Interpolate() gives a smoother evaluation than raw bin lookup, useful
// since the fit x-values won't generally sit exactly on bin centers.
double PSShape(double x)
{
    double xmin = g_hPS->GetXaxis()->GetXmin();
    double xmax = g_hPS->GetXaxis()->GetXmax();
    if(x < xmin || x > xmax)
        return 0.0;
    return g_hPS->Interpolate(x);
}

// Total fit function: Amplitude_sig * convolution(as,r0,eps,sigma_res) + Amplitude_PS * PS(x)
// par: [0]=Amp_sig, [1]=as, [2]=r0, [3]=eps, [4]=sigma_res, [5]=Amp_PS
double FitFunction(double* x, double* par)
{
    double amp_sig = par[0];
    double a_s = par[1];
    double r0 = par[2];
    double eps = par[3];
    double sigma_res = par[4];
    double amp_ps = par[5];

    double p_conv[5] = {a_s, r0, eps, 0.0, sigma_res}; // {as, r0, eps} + {mean=0, sigma}
    g_f_conv->SetParameters(p_conv);
    double sig = amp_sig * g_f_conv->Eval(x[0]);

    double ps = amp_ps * PSShape(x[0]);

    return sig + ps;
}

void fit_spectra_12Li_sWave()
{
    double fitmin = -1.0;
    double fitmax = 8.0;

    // Ex spectra histogram file and name
    std::string filename = "./Inputs/Ex_total_heavyGate_11Li_d_p.root";
    // const char* histname = "hEx_9Li";
    std::string histname = "hEx_Total";

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);

    // --- Read the data histogram ---
    TFile* file = TFile::Open(filename.c_str(), "READ");
    if(!file || file->IsZombie())
    {
        std::cerr << "Could not open file: " << filename << std::endl;
        return;
    }

    TH1* h = dynamic_cast<TH1*>(file->Get(histname.c_str()));
    if(!h)
    {
        std::cerr << "Histogram '" << histname << "' was not found in " << filename << std::endl;
        return;
    }
    h->SetDirectory(nullptr); // Detach from the file so it can be closed
    file->Close();

    // Phase space
    ROOT::RDataFrame phase {"SimulationTTree", "../../Simulation/Outputs/7Li/2H_1H_TRIUMF_Eex_0.000_nPS_1_pPS_0.root"};
    // Copy exactly the same binning as the experimental histogram
    auto model =
        ROOT::RDF::TH1DModel("hPS", "Phase Space", h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());

    auto hPS = phase.Histo1D(model, "EexGateHeavy_side", "weight");

    g_hPS = (TH1D*)hPS.GetPtr()->Clone("hPSClone");
    g_hPS->SetDirectory(nullptr);

    //--------------------------------------------------------------
    // Normalize
    //--------------------------------------------------------------
    double integral = g_hPS->Integral();

    if(integral > 0)
        g_hPS->Scale(1.0 / integral);
    else
        std::cerr << "PS histogram has zero integral!" << std::endl;

    //--------------------------------------------------------------
    // Draw only the PS histogram
    //--------------------------------------------------------------
    auto* c2 = new TCanvas("cPS", "Phase Space", 800, 600);

    g_hPS->SetLineColor(kBlue + 1);
    g_hPS->SetLineWidth(2);
    g_hPS->SetFillStyle(0);
    g_hPS->SetTitle("Phase Space;E_{nf} (MeV);Arbitrary units");

    hPS->DrawClone("hist");

    // --- Build the s-wave convolution (wider internal range than the fit
    //     range, to avoid truncating the Gaussian resolution tails) ---
    const double conv_min = fitmin - 0.6;
    const double conv_max = fitmax + 0.6;

    g_f_orig = new TF1("g_f_orig", SWaveSpectrum, conv_min, conv_max, 3);
    g_f_gaus = new TF1("g_f_gaus", "TMath::Gaus(x,[0],[1],true)", conv_min, conv_max);
    g_conv = new TF1Convolution(g_f_orig, g_f_gaus, conv_min, conv_max, true);
    g_conv->SetNofPointsFFT(1000);
    g_f_conv = new TF1("g_f_conv", *g_conv, conv_min, conv_max, g_conv->GetNpar());

    // --- Fit function ---
    TF1* fit_func = new TF1("fit_func", FitFunction, fitmin, fitmax, 6);
    fit_func->SetParNames("Amplitude_{sig}", "a_{s} (fm)", "r0 (fm)", "#varepsilon (MeV)", "#sigma_{res} (MeV)",
                          "Amplitude_{PS}");

    // Initial values -- adjust the amplitude guesses to the scale of your histogram
    double amp_guess = h->GetMaximum() > 0 ? h->GetMaximum() : 1.0;
    fit_func->SetParameter(0, amp_guess);       // Amplitude_sig, free
    fit_func->SetParameter(1, -20.0);           // as, free (example initial value)
    fit_func->SetParameter(2, 3.0);             // r0, fixed later
    fit_func->SetParameter(3, 0.370);           // eps, fixed later
    fit_func->SetParameter(4, 0.24);            // sigma_res, fixed later
    fit_func->SetParameter(5, 0.5 * amp_guess); // Amplitude_PS, free (rough starting guess)

    // Reasonable limits for the free parameters (adjust if needed)
    fit_func->SetParLimits(0, 0.0, 1e6 * amp_guess);
    fit_func->SetParLimits(1, -500, -0.5); // Scattering length is typically negative
    fit_func->FixParameter(5, 0.0);

    // Fix all parameters except Amplitude_sig, as, and Amplitude_PS
    fit_func->FixParameter(2, 3.0);   // r0
    fit_func->FixParameter(3, 0.370); // eps (approximately the 2n separation energy of 11Li)
    fit_func->FixParameter(4, 0.24);  // sigma_res

    // --- Perform the fit ---
    // "R" = use the function range, "S" = return a TFitResultPtr.
    // Add "L" if the histogram has low statistics (likelihood fit
    // instead of chi-square), e.g. "RSL".
    TFitResultPtr fit_result = h->Fit(fit_func, "RS");

    std::cout << "\n=== Fit results ===\n";
    std::cout << "a_s = " << fit_func->GetParameter(1) << " +/- " << fit_func->GetParError(1) << " fm\n";
    std::cout << "Amplitude_sig = " << fit_func->GetParameter(0) << " +/- " << fit_func->GetParError(0) << "\n";
    std::cout << "Amplitude_PS  = " << fit_func->GetParameter(5) << " +/- " << fit_func->GetParError(5) << "\n";
    std::cout << "chi2/ndf = " << fit_func->GetChisquare() << " / " << fit_func->GetNDF() << "\n";

    // --- Draw the results: data, total fit, and the two components separately ---
    TCanvas* c1 = new TCanvas("c", "Fit s-wave + PS", 800, 600);
    h->SetTitle("s-wave + Phase Space fit;E_{nf} (MeV);Counts");
    h->Draw("E");

    fit_func->SetLineColor(kRed + 1);
    fit_func->SetLineWidth(2);
    fit_func->Draw("SAME");

    // Signal-only component, using the fitted parameters
    TF1* f_sig_only = new TF1(
        "f_sig_only",
        [](double* xx, double* pp)
        {
            double p_conv[5] = {pp[1], pp[2], pp[3], 0.0, pp[4]};
            g_f_conv->SetParameters(p_conv);
            return pp[0] * g_f_conv->Eval(xx[0]);
        },
        fitmin, fitmax, 6);
    f_sig_only->SetParameters(fit_func->GetParameters());
    f_sig_only->SetLineColor(kBlue + 1);
    f_sig_only->SetLineStyle(2);
    f_sig_only->SetLineWidth(2);
    f_sig_only->Draw("SAME");

    // PS-only component, using the fitted amplitude
    TF1* f_ps_only =
        new TF1("f_ps_only", [](double* xx, double* pp) { return pp[5] * PSShape(xx[0]); }, fitmin, fitmax, 6);
    f_ps_only->SetParameters(fit_func->GetParameters());
    f_ps_only->SetLineColor(kGreen + 2);
    f_ps_only->SetLineStyle(2);
    f_ps_only->SetLineWidth(2);
    f_ps_only->Draw("SAME");

    TLegend* leg = new TLegend(0.55, 0.65, 0.88, 0.88);
    leg->AddEntry(h, "Data", "lep");
    leg->AddEntry(fit_func, "Total fit", "l");
    leg->AddEntry(f_sig_only, "s-wave component", "l");
    leg->AddEntry(f_ps_only, "Phase Space component", "l");
    leg->Draw();

    c1->SaveAs("fit_s_wave_histogram.png");
    c1->SaveAs("fit_s_wave_histogram.pdf");
}