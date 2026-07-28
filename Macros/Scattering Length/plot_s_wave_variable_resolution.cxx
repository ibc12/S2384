// ---------------------------------------------------------------------
// Comparacion: convolucion con resolucion FIJA (sigma=0.14 MeV, via
// TF1Convolution) vs "smearing" con resolucion VARIABLE en energia
// (sigma(E) interpolada linealmente entre sigma(0)=0.14 y sigma(1)=0.16 MeV),
// calculada mediante integracion numerica punto a punto.
//
// IMPORTANTE: una resolucion que depende de la energia no es una
// convolucion en sentido estricto (falta invariancia traslacional), por
// eso TF1Convolution no sirve para el caso variable. Aqui se calcula
// "a mano" con una suma tipo trapecio.
//
// Convencion elegida: sigma se evalua en la energia VERDADERA E' (la
// variable de integracion), es decir cada evento fisico se difumina con
// la resolucion propia de su energia. Si prefieres evaluar sigma en la
// energia observada Eobs, cambia SigmaRes(Ep) -> SigmaRes(Eobs) en
// VariableSigmaSmear (dejo un comentario en ese punto).
//
// Uso en ROOT:
//   root -l plot_s_wave_variable_resolution.C
// ---------------------------------------------------------------------

#include <TF1.h>
#include <TF1Convolution.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TMath.h>

static const double HBAR_C = 197.3269804; // MeV*fm

double MuC2_8He_n() {
    const double m_n_u    = 1.008665;
    const double m_8He_u  = 8.033934;
    const double u_to_MeV = 931.494;
    return (m_n_u * m_8He_u) / (m_n_u + m_8He_u) * u_to_MeV;
}
static const double MU_C2 = MuC2_8He_n();

double WaveNumber(double E) {
    return TMath::Sqrt(2.0 * MU_C2 * E) / HBAR_C;
}

// Funcion fisica (Ecs. 1-2), 0 para Enf<=0. Parametros: [0]=as,[1]=r0,[2]=eps
double SWaveSpectrum(double *x, double *par) {
    double Enf = x[0];
    if (Enf <= 0.0) return 0.0;

    double a_s = par[0];
    double r0  = par[1];
    double eps = par[2];

    double knf = WaveNumber(Enf);
    if (knf < 1e-8) knf = 1e-8;

    double k = WaveNumber(eps);

    double cot_delta = (-1.0 / a_s + 0.5 * r0 * knf * knf) / knf;
    double sin_delta = 1.0 / TMath::Sqrt(1.0 + cot_delta * cot_delta);
    double cos_delta = cot_delta / TMath::Sqrt(1.0 + cot_delta * cot_delta);

    double bracket = (1.0 / (k * k + knf * knf)) * cos_delta
                    + (k / knf) * sin_delta;

    return knf * bracket * bracket;
}

// Resolucion dependiente de la energia: interpolacion lineal entre
// sigma(0 MeV)=0.14 y sigma(1 MeV)=0.16, extrapolada linealmente sin
// limite fuera de ese rango (misma pendiente en todo el dominio).
double SigmaRes(double E) {
    const double E0 = 0.0, s0 = 0.14;
    const double E1 = 1.0, s1 = 0.2;
    return s0 + (s1 - s0) * (E - E0) / (E1 - E0);
}

// Integral de convolucion con sigma variable, calculada con regla del
// trapecio sobre una malla fina de la energia verdadera E'.
double VariableSigmaSmear(double *x, double *par) {
    double a_s = par[0];
    double r0  = par[1];
    double eps = par[2];
    double Eobs = x[0];

    const double Ep_min = -0.3, Ep_max = 2.5; // rango de integracion (energia verdadera)
    const int N = 800;
    const double dE = (Ep_max - Ep_min) / N;

    double parPhys[3] = {a_s, r0, eps};
    double sum = 0.0;

    for (int j = 0; j <= N; ++j) {
        double Ep = Ep_min + j * dE;
        double xx[1] = {Ep};
        double fphys = SWaveSpectrum(xx, parPhys);
        if (fphys == 0.0) continue;

        double sigma = SigmaRes(Ep);  // <-- convencion: sigma en la energia VERDADERA
        double kernel = TMath::Gaus(Eobs - Ep, 0.0, sigma, true);

        double weight = (j == 0 || j == N) ? 0.5 : 1.0;
        sum += weight * fphys * kernel;
    }
    return sum * dE;
}

void plot_s_wave_variable_resolution() {
    gStyle->SetOptStat(0);

    const double r0  = 3.0;
    const double eps = 0.20;
    const double sigma_fixed = 0.14; // MeV, resolucion fija de referencia

    const double xmin = -0.6, xmax = 2.5;   // rango interno (para TF1Convolution)
    const double drawmin = -0.3, drawmax = 2.0;

    double as_values[3] = {-10.0, -20.0, -50.0};
    Color_t colors[3] = {kBlue + 1, kRed + 1, kGreen + 2};

    TCanvas *c1 = new TCanvas("c1", "Resolucion fija vs variable", 1500, 500);
    c1->Divide(3, 1);

    static TF1 *f_orig[3];
    static TF1 *f_gaus[3];
    static TF1Convolution *conv[3];
    static TF1 *f_conv_fixed_raw[3];
    static TF1 *f_orig_n[3];
    static TF1 *f_conv_fixed_n[3];
    static TF1 *f_conv_var_n[3];
    static TLegend *leg[3];

    for (int i = 0; i < 3; ++i) {
        double a_s = as_values[i];

        // --- funcion fisica original ---
        f_orig[i] = new TF1(Form("f_orig_%d", i), SWaveSpectrum, xmin, xmax, 3);
        f_orig[i]->SetParameters(a_s, r0, eps);
        f_orig[i]->SetNpx(2000);

        // --- convolucion con sigma FIJA (TF1Convolution, como antes) ---
        f_gaus[i] = new TF1(Form("f_gaus_%d", i), "TMath::Gaus(x,[0],[1],true)", xmin, xmax);
        f_gaus[i]->SetParameters(0.0, sigma_fixed);

        conv[i] = new TF1Convolution(f_orig[i], f_gaus[i], xmin, xmax, true);
        conv[i]->SetNofPointsFFT(1000);

        f_conv_fixed_raw[i] = new TF1(Form("f_conv_fixed_raw_%d", i), *conv[i], xmin, xmax, conv[i]->GetNpar());
        f_conv_fixed_raw[i]->SetParameters(a_s, r0, eps, 0.0, sigma_fixed);
        f_conv_fixed_raw[i]->SetNpx(2000);

        // --- smearing con sigma VARIABLE (integracion numerica) ---
        TF1 *f_conv_var_raw = new TF1(Form("f_conv_var_raw_%d", i), VariableSigmaSmear, xmin, xmax, 3);
        f_conv_var_raw->SetParameters(a_s, r0, eps);
        f_conv_var_raw->SetNpx(500); // mas caro de evaluar: menos puntos

        // --- normalizar las tres curvas al maximo de la original ---
        double ymax = f_orig[i]->GetMaximum(0.0, drawmax);

        f_orig_n[i] = new TF1(Form("f_orig_n_%d", i),
            [=](double *xx, double *pp) { return SWaveSpectrum(xx, pp) / ymax; },
            drawmin, drawmax, 3);
        f_orig_n[i]->SetParameters(a_s, r0, eps);
        f_orig_n[i]->SetNpx(2000);
        f_orig_n[i]->SetLineColor(colors[i]);
        f_orig_n[i]->SetLineWidth(2);
        f_orig_n[i]->SetLineStyle(1);

        TF1 *f_conv_fixed_i = f_conv_fixed_raw[i];
        f_conv_fixed_n[i] = new TF1(Form("f_conv_fixed_n_%d", i),
            [=](double *xx, double *pp) { return f_conv_fixed_i->Eval(xx[0]) / ymax; },
            drawmin, drawmax, 0);
        f_conv_fixed_n[i]->SetNpx(2000);
        f_conv_fixed_n[i]->SetLineColor(colors[i]);
        f_conv_fixed_n[i]->SetLineWidth(2);
        f_conv_fixed_n[i]->SetLineStyle(2);

        f_conv_var_n[i] = new TF1(Form("f_conv_var_n_%d", i),
            [=](double *xx, double *pp) { return VariableSigmaSmear(xx, pp) / ymax; },
            drawmin, drawmax, 3);
        f_conv_var_n[i]->SetParameters(a_s, r0, eps);
        f_conv_var_n[i]->SetNpx(400);
        f_conv_var_n[i]->SetLineColor(colors[i]);
        f_conv_var_n[i]->SetLineWidth(2);
        f_conv_var_n[i]->SetLineStyle(3);

        c1->cd(i + 1);
        f_orig_n[i]->SetTitle(Form("a_{s} = %.0f fm;E_{nf} (MeV);d#sigma/dE_{nf} (u. arb.)", a_s));
        f_orig_n[i]->SetMinimum(0.0);
        f_orig_n[i]->SetMaximum(1.15);
        f_orig_n[i]->Draw();
        f_conv_fixed_n[i]->Draw("SAME");
        f_conv_var_n[i]->Draw("SAME");

        leg[i] = new TLegend(0.35, 0.65, 0.88, 0.88);
        leg[i]->AddEntry(f_orig_n[i], "Sin convolucionar", "l");
        leg[i]->AddEntry(f_conv_fixed_n[i], Form("#sigma fija = %.2f MeV", sigma_fixed), "l");
        leg[i]->AddEntry(f_conv_var_n[i], "#sigma(E)=0.14#rightarrow0.2 MeV", "l");
        leg[i]->Draw();
    }

    c1->SaveAs("s_wave_variable_resolution.png");
    c1->SaveAs("s_wave_variable_resolution.pdf");
}