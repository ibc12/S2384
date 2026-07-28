// ---------------------------------------------------------------------
// Espectro de energia relativa dsigma/dEnf para dispersion s-wave (l=0)
// en el sistema 8He + n, segun las Ecs. (1) y (2) de
// H.T. Johansson et al., Nucl. Phys. A 842 (2010) 15-32.
//
// Ec. (2):  knf * cot(delta) = -1/as + (1/2) r0 knf^2   (se ignora O(knf^4))
// Ec. (1):  dsigma/dEnf ~ knf * [ 1/(k^2+knf^2)*cos(delta) + (k/knf)*sin(delta) ]^2
//
// knf = sqrt(2*mu*Enf)/hbar   -> numero de onda asociado a la energia del eje x
// k   = sqrt(2*mu*eps)/hbar   -> parametro fijo definido por epsilon (~ S2n en el fit)
//
// Uso en ROOT:
//   root -l plot_s_wave_resonance.C
// ---------------------------------------------------------------------

#include <TAxis.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLegend.h>
#include <TMath.h>
#include <TStyle.h>

static const double HBAR_C = 197.3269804; // MeV*fm

// masa reducida 8He+n en MeV/c^2
double ReducedMass_8He_n()
{
    const double m_n_u = 1.008665;
    const double m_8He_u = 8.033934;
    const double u_to_MeV = 931.494;
    return (m_n_u * m_8He_u) / (m_n_u + m_8He_u) * u_to_MeV;
}

// numero de onda k(E) = sqrt(2*mu*E)/hbar  [fm^-1], E en MeV
double WaveNumber(double E, double mu_c2)
{
    return TMath::Sqrt(2.0 * mu_c2 * E) / HBAR_C;
}

// Funcion del espectro para usar con TF1.
// Parametros: [0]=as (fm), [1]=r0 (fm), [2]=eps (MeV), [3]=mu_c2 (MeV)
double SWaveSpectrum(double* x, double* par)
{
    double Enf = x[0];
    if(Enf <= 0.0)
        return 0.0;

    double a_s = par[0];
    double r0 = par[1];
    double eps = par[2];
    double mu_c2 = par[3];

    double knf = WaveNumber(Enf, mu_c2);
    if(knf < 1e-8)
        knf = 1e-8; // evita division por cero cerca de Enf=0

    double k = WaveNumber(eps, mu_c2);

    // Ec. (2), sin el termino O(knf^4)
    double cot_delta = (-1.0 / a_s + 0.5 * r0 * knf * knf) / knf;

    double sin_delta = 1.0 / TMath::Sqrt(1.0 + cot_delta * cot_delta);
    double cos_delta = cot_delta / TMath::Sqrt(1.0 + cot_delta * cot_delta);

    double bracket = (1.0 / (k * k + knf * knf)) * cos_delta + (k / knf) * sin_delta;

    return knf * bracket * bracket;
}

void plot_s_wave_resonance()
{
    gStyle->SetOptStat(0);

    double mu_c2 = ReducedMass_8He_n();

    // Parametros de ejemplo (editables)
    double r0 = 3.0;   // fm, fijo segun el paper
    double eps = 0.20; // MeV

    const double Emin = -1e-4, Emax = 2.0; // MeV

    TCanvas* c1 = new TCanvas("c1", "Espectro s-wave 8He+n", 800, 600);

    // Curva principal: as = -20 fm
    TF1* f_main = new TF1("f_main", SWaveSpectrum, Emin, Emax, 4);
    f_main->SetParameters(-20.0, r0, eps, mu_c2);
    f_main->SetNpx(2000);
    f_main->SetLineColor(kBlue + 1);
    f_main->SetLineWidth(2);

    // Curvas de comparacion variando as
    TF1* f_alt1 = new TF1("f_alt1", SWaveSpectrum, Emin, Emax, 4);
    f_alt1->SetParameters(-10.0, r0, eps, mu_c2);
    f_alt1->SetNpx(2000);
    f_alt1->SetLineColor(kRed + 1);
    f_alt1->SetLineStyle(2);
    f_alt1->SetLineWidth(2);

    TF1* f_alt2 = new TF1("f_alt2", SWaveSpectrum, Emin, Emax, 4);
    f_alt2->SetParameters(-50.0, r0, eps, mu_c2);
    f_alt2->SetNpx(2000);
    f_alt2->SetLineColor(kGreen + 2);
    f_alt2->SetLineStyle(2);
    f_alt2->SetLineWidth(2);

    // Normalizar cada curva a su maximo evaluando en una malla fina
    auto normalize = [&](TF1* f)
    {
        double ymax = 0.0;
        for(double E = Emin; E < Emax; E += (Emax - Emin) / 2000.0)
        {
            double y = f->Eval(E);
            if(y > ymax)
                ymax = y;
        }
        // Guardamos el factor de normalizacion como parametro extra via escala manual
        return ymax;
    };

    double ymax_main = normalize(f_main);
    double ymax_alt1 = normalize(f_alt1);
    double ymax_alt2 = normalize(f_alt2);

    // Como TF1 no soporta division directa por una constante calculada,
    // usamos TF1 tipo "NSUM"/wrapper con una escala explicita.
    TF1* f_main_n =
        new TF1("f_main_n", [=](double* x, double* par) { return SWaveSpectrum(x, par) / ymax_main; }, Emin, Emax, 4);
    f_main_n->SetParameters(-20.0, r0, eps, mu_c2);
    f_main_n->SetNpx(2000);
    f_main_n->SetLineColor(kBlue + 1);
    f_main_n->SetLineWidth(2);

    TF1* f_alt1_n =
        new TF1("f_alt1_n", [=](double* x, double* par) { return SWaveSpectrum(x, par) / ymax_alt1; }, Emin, Emax, 4);
    f_alt1_n->SetParameters(-10.0, r0, eps, mu_c2);
    f_alt1_n->SetNpx(2000);
    f_alt1_n->SetLineColor(kRed + 1);
    f_alt1_n->SetLineStyle(2);
    f_alt1_n->SetLineWidth(2);

    TF1* f_alt2_n =
        new TF1("f_alt2_n", [=](double* x, double* par) { return SWaveSpectrum(x, par) / ymax_alt2; }, Emin, Emax, 4);
    f_alt2_n->SetParameters(-50.0, r0, eps, mu_c2);
    f_alt2_n->SetNpx(2000);
    f_alt2_n->SetLineColor(kGreen + 2);
    f_alt2_n->SetLineStyle(2);
    f_alt2_n->SetLineWidth(2);

    f_main_n->SetTitle(
        "Espectro de energia relativa 8He+n (onda s, Ecs. 1-2);E_{nf} (MeV);d#sigma/dE_{nf} (u. arb., normalizado)");
    f_main_n->SetMinimum(0.0);
    f_main_n->SetMaximum(1.15);
    f_main_n->Draw();
    f_alt1_n->Draw("SAME");
    f_alt2_n->Draw("SAME");

    TLegend* leg = new TLegend(0.55, 0.65, 0.88, 0.88);
    leg->AddEntry(f_main_n, Form("a_{s}=-20 fm, r0=%.1f fm, #varepsilon=%.2f MeV", r0, eps), "l");
    leg->AddEntry(f_alt1_n, Form("a_{s}=-10 fm, r0=%.1f fm, #varepsilon=%.2f MeV", r0, eps), "l");
    leg->AddEntry(f_alt2_n, Form("a_{s}=-50 fm, r0=%.1f fm, #varepsilon=%.2f MeV", r0, eps), "l");
    leg->Draw();

    c1->SaveAs("s_wave_spectrum_ROOT.png");
    c1->SaveAs("s_wave_spectrum_ROOT.pdf");
}
