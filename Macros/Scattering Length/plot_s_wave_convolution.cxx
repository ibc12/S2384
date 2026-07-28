// ---------------------------------------------------------------------
// Convolucion del espectro de energia relativa s-wave (8He+n) con la
// resolucion experimental, modelada como una gaussiana de sigma = 0.14 MeV.
//
// Ecs. (1) y (2) de H.T. Johansson et al., Nucl. Phys. A 842 (2010) 15-32.
//
// Se usa TF1Convolution para convolucionar la funcion fisica (SWaveSpectrum)
// con el kernel gaussiano de resolucion. Canvas dividido en 3 pads, uno por
// cada valor de a_s, mostrando la curva original y su convolucion.
//
// Uso en ROOT:
//   root -l plot_s_wave_convolution.C
// ---------------------------------------------------------------------

#include <TCanvas.h>
#include <TF1.h>
#include <TF1Convolution.h>
#include <TLegend.h>
#include <TMath.h>
#include <TStyle.h>

static const double HBAR_C = 197.3269804; // MeV*fm

// Masa reducida 8He+n en MeV/c^2 (constante global, no es parametro libre)
double MuC2_8He_n()
{
    const double m_n_u = 1.008665;
    const double m_8He_u = 8.033934;
    const double u_to_MeV = 931.494;
    return (m_n_u * m_8He_u) / (m_n_u + m_8He_u) * u_to_MeV;
}
static const double MU_C2 = MuC2_8He_n();

// numero de onda k(E) = sqrt(2*mu*E)/hbar  [fm^-1], E en MeV
// se usa abs(E) para permitir un rango ligeramente negativo requerido
// por la cola de la gaussiana en la convolucion (E<0 no es fisico, pero
// mantiene la funcion continua/estable cerca del umbral E=0)
double WaveNumber(double E)
{
    return TMath::Sqrt(2.0 * MU_C2 * E) / HBAR_C;
}

// Funcion fisica del espectro (Ecs. 1-2), sin el termino O(knf^4).
// Parametros: [0]=as (fm), [1]=r0 (fm), [2]=eps (MeV)
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
        knf = 1e-8; // evita division por cero cerca de Enf=0

    double k = WaveNumber(eps);

    double cot_delta = (-1.0 / a_s + 0.5 * r0 * knf * knf) / knf;
    double sin_delta = 1.0 / TMath::Sqrt(1.0 + cot_delta * cot_delta);
    double cos_delta = cot_delta / TMath::Sqrt(1.0 + cot_delta * cot_delta);

    double bracket = (1.0 / (k * k + knf * knf)) * cos_delta + (k / knf) * sin_delta;

    return knf * bracket * bracket;
}

void plot_s_wave_convolution()
{
    gStyle->SetOptStat(0);

    const double r0 = 3.0;         // fm, fijo segun el paper
    const double eps = 0.20;       // MeV
    const double sigma_res = 0.14; // MeV, resolucion experimental

    // Rango de la funcion / convolucion. Se extiende un poco por debajo de
    // 0 para no cortar la cola izquierda de la gaussiana de resolucion.
    const double xmin = -0.6, xmax = 2.5;
    const double drawmin = -0.6, drawmax = 2.0; // rango fisico a dibujar

    double as_values[3] = {-10.0, -20.0, -50.0};
    Color_t colors[3] = {kBlue + 1, kRed + 1, kGreen + 2};

    TCanvas* c1 = new TCanvas("c1", "s-wave y convolucion con resolucion", 1500, 500);
    c1->Divide(3, 1);

    // Punteros que deben sobrevivir fuera del bucle (ROOT no copia TF1 al dibujar)
    static TF1* f_orig[3];
    static TF1* f_gaus[3];
    static TF1Convolution* conv[3];
    static TF1* f_conv[3];
    static TLegend* leg[3];

    for(int i = 0; i < 3; ++i)
    {
        double a_s = as_values[i];

        // --- funcion original (sin convolucionar) ---
        f_orig[i] = new TF1(Form("f_orig_%d", i), SWaveSpectrum, xmin, xmax, 3);
        f_orig[i]->SetParameters(a_s, r0, eps);
        f_orig[i]->SetNpx(2000);

        // --- kernel gaussiano de resolucion, normalizado a integral 1 ---
        // par[0] = media (fija en 0), par[1] = sigma
        f_gaus[i] = new TF1(Form("f_gaus_%d", i), "TMath::Gaus(x,[0],[1],true)", xmin, xmax);
        f_gaus[i]->SetParameters(0.0, sigma_res);

        // --- convolucion s-wave (x) gaussiana ---
        conv[i] = new TF1Convolution(f_orig[i], f_gaus[i], xmin, xmax, true);
        conv[i]->SetNofPointsFFT(1000);

        f_conv[i] = new TF1(Form("f_conv_%d", i), *conv[i], xmin, xmax, conv[i]->GetNpar());
        // parametros: primero los de f_orig (as, r0, eps), luego los de f_gaus (media, sigma)
        f_conv[i]->SetParameters(a_s, r0, eps, 0.0, sigma_res);
        f_conv[i]->SetNpx(2000);

        // --- normalizar ambas curvas al maximo de la funcion original ---
        double ymax = f_orig[i]->GetMaximum(drawmin, drawmax);

        TF1* f_orig_n = new TF1(
            Form("f_orig_n_%d", i), [=](double* xx, double* pp) { return SWaveSpectrum(xx, pp) / ymax; }, drawmin,
            drawmax, 3);
        f_orig_n->SetParameters(a_s, r0, eps);
        f_orig_n->SetNpx(2000);
        f_orig_n->SetLineColor(colors[i]);
        f_orig_n->SetLineWidth(2);
        f_orig_n->SetLineStyle(1);

        TF1* f_conv_n = new TF1(
            Form("f_conv_n_%d", i), [=](double* xx, double* pp) { return f_conv[i]->Eval(xx[0]) / ymax; }, drawmin,
            drawmax, 0);
        f_conv_n->SetNpx(2000);
        f_conv_n->SetLineColor(colors[i]);
        f_conv_n->SetLineWidth(2);
        f_conv_n->SetLineStyle(2);

        c1->cd(i + 1);
        f_orig_n->SetTitle(Form("a_{s} = %.0f fm;E_{nf} (MeV);d#sigma/dE_{nf} (u. arb.)", a_s));
        f_orig_n->SetMinimum(0.0);
        f_orig_n->SetMaximum(1.15);
        f_orig_n->Draw();
        f_conv_n->Draw("SAME");

        leg[i] = new TLegend(0.45, 0.72, 0.88, 0.88);
        leg[i]->AddEntry(f_orig_n, "Sin convolucionar", "l");
        leg[i]->AddEntry(f_conv_n, Form("Convol. #sigma_{res}=%.2f MeV", sigma_res), "l");
        leg[i]->Draw();
    }

    c1->SaveAs("s_wave_convolution.png");
    c1->SaveAs("s_wave_convolution.pdf");
}