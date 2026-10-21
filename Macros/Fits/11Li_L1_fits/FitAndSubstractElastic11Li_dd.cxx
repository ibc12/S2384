// Macro: ajuste del estado fundamental en el espectro de Ex y sustraccion de su contribucion
//
// Uso:   root -l fit_gs_Ex.C
//   (o)  root -l -b -q fit_gs_Ex.C
//
// Nota: se asume que ActRoot::MergerData tiene el miembro fEx y que la
// columna del TTree se llama "MergerData". Ajustar si fuera necesario.

#include "ActMergerData.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TFile.h"

#include <cmath>
#include <iostream>

#include "../../../PrettyStyle.C"

void FitAndSubstractElastic11Li_dd()
{
    PrettyStyle(true, false);
    // ---------------- Configuracion ----------------
    const std::string fileName {"../../../PostAnalysis/Outputs/tree_ex_F_11Li_d_d_filtered.root"};
    const std::string treeName {"Final_Tree"};

    const double binWidth {0.150}; // ancho de bin en MeV (150 keV)
    const double xMin {-3.};       // MeV
    const double xMax {12.};       // MeV (rango / binWidth debe ser entero: 15 / 0.15 = 100 bins)
    const int nBins {static_cast<int>(std::lround((xMax - xMin) / binWidth))};

    // Titulo del eje Y con el ancho de bin en keV
    const std::string yTitle {Form("Counts / %.0f keV", binWidth * 1000.)};

    const double fitMin {-1.};  // inicio del rango de ajuste (MeV)
    const double fitMax {0.6};  // fin del rango de ajuste (MeV): Ex < 0.3
    const bool fixMean {false}; // true -> fija la media en Ex = 0
    // -------------------------------------------------

    gStyle->SetOptStat(0);

    // ---------------- RDataFrame -------------------
    ROOT::RDataFrame df {treeName, fileName};
    auto def {df.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"})};

    const std::string hTitle {"Espectro de Ex;E_{x} [MeV];" + yTitle};
    auto hExPtr {def.Histo1D({"hEx", hTitle.c_str(), nBins, xMin, xMax}, "Ex")};
    TH1D* hEx {(TH1D*)hExPtr->Clone("hEx_total")};
    hEx->SetDirectory(nullptr);

    // ---------------- Ajuste gaussiano -------------
    // Valores iniciales: altura = maximo en la region del fit, media = 0, sigma = 0.15 MeV
    hEx->GetXaxis()->SetRangeUser(fitMin, fitMax);
    double amp0 {hEx->GetMaximum()};
    hEx->GetXaxis()->SetRange(0, 0); // resetea el rango

    auto* gs {new TF1("gs", "gaus", fitMin, fitMax)};
    gs->SetParNames("Amp", "Mean", "Sigma");
    gs->SetParameters(amp0, 0., 0.15);
    if(fixMean)
        gs->FixParameter(1, 0.);
    else
        gs->SetParLimits(1, -0.3, 0.3);
    gs->SetParLimits(2, 0.01, 1.);
    gs->SetLineColor(kRed);
    gs->SetLineWidth(2);

    // R: usa el rango de la funcion | L: likelihood (Poisson, bueno con pocas cuentas)
    // I: integra la funcion dentro de cada bin | S: devuelve el resultado
    TFitResultPtr res {hEx->Fit(gs, "RLIS")};

    const double amp {gs->GetParameter(0)};
    const double mean {gs->GetParameter(1)};
    const double sigma {gs->GetParameter(2)};
    const double fwhm {2.3548 * sigma};
    const double binW {hEx->GetBinWidth(1)};
    const double nGS {amp * sigma * std::sqrt(2. * M_PI) / binW}; // cuentas totales de la gaussiana

    std::cout << "---- Ajuste del estado fundamental ----" << '\n';
    std::cout << "Media  = " << mean << " MeV" << '\n';
    std::cout << "Sigma  = " << sigma << " MeV  (FWHM = " << fwhm << " MeV)" << '\n';
    std::cout << "Cuentas del gs (integral gaussiana) = " << nGS << '\n';

    // ---------------- Sustraccion bin a bin --------
    // Para cada bin se resta el valor esperado de la gaussiana en ese bin,
    // calculado como integral de la funcion en el bin / ancho del bin.
    // Se usa una gaussiana extendida a todo el rango (no solo al rango del fit).
    auto* gsFull {new TF1("gsFull", "gaus", xMin, xMax)};
    gsFull->SetParameters(amp, mean, sigma);

    auto* hSub {(TH1D*)hEx->Clone("hEx_sub")};
    hSub->SetTitle(("Espectro de Ex sin estado fundamental;E_{x} [MeV];" + yTitle).c_str());
    hSub->SetDirectory(nullptr);
    hSub->GetListOfFunctions()->Clear(); // quita el fit heredado del clon

    for(int i = 1; i <= hSub->GetNbinsX(); i++)
    {
        const double lo {hSub->GetXaxis()->GetBinLowEdge(i)};
        const double up {hSub->GetXaxis()->GetBinUpEdge(i)};
        const double expected {gsFull->Integral(lo, up) / hSub->GetBinWidth(i)};

        const double content {hEx->GetBinContent(i)};
        const double err {hEx->GetBinError(i)};

        hSub->SetBinContent(i, content - expected);
        // error: estadistico del bin original + fluctuacion poissoniana de lo restado
        hSub->SetBinError(i, std::sqrt(err * err + expected));
        // Si prefieres evitar bins negativos: hSub->SetBinContent(i, std::max(0., content - expected));
    }

    // ---------------- Dibujo -----------------------
    auto* c {new TCanvas("cEx", "Ex: total vs sin gs", 1100, 700)};
    c->SetGrid();

    hEx->SetLineColor(kBlack);
    hEx->SetLineWidth(2);
    hEx->SetTitle(("Espectro de E_{x};E_{x} [MeV];" + yTitle).c_str());

    hSub->SetLineColor(kAzure + 1);
    hSub->SetLineWidth(2);
    hSub->SetFillColorAlpha(kAzure + 1, 0.25);

    hEx->Draw("HIST E");
    hSub->Draw("HIST E SAME");
    gsFull->SetLineColor(kRed);
    gsFull->SetLineStyle(2);
    gsFull->SetLineWidth(2);
    gsFull->Draw("SAME");

    auto* leg {new TLegend(0.55, 0.62, 0.89, 0.88)};
    leg->SetBorderSize(0);
    leg->AddEntry(hEx, "Total", "l");
    leg->AddEntry(gsFull, Form("Fit gauss. gs (#sigma = %.0f keV)", sigma * 1000.), "l");
    leg->AddEntry(hSub, "Total - gs", "lf");
    leg->Draw();

    auto* c2 {new TCanvas("cExSub", "Ex sin estado fundamental", 1100, 700)};
    c2->SetGrid();

    hSub->SetTitle(("Espectro de E_{x} sin estado fundamental;E_{x} [MeV];" + yTitle).c_str());
    hSub->SetFillStyle(0);                     // sin relleno
    hSub->GetXaxis()->SetRangeUser(-1, 10); // opcional
    hSub->GetYaxis()->SetRangeUser(0, 1.2 * hSub->GetMaximum()); // opcional
    hSub->Draw("HIST E");

    c2->Modified();
    c2->Update();

    auto* fOut {TFile::Open("./Outputs/Ex_sin_gs.root", "RECREATE")};
    fOut->cd();
    hSub->Write("hEx_sub");
    hEx->Write("hEx_total");
    fOut->Close();
}