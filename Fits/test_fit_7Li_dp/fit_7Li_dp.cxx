#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TSpline.h"
#include "TMath.h"

#include "ROOT/RDataFrame.hxx"

TSpline3 *gSplinePS = nullptr;

// función que evalúa la spline + gaussianas
double fitFunc(double *x, double *p) {
    double xx = x[0];

    // --- gaussianas existentes ---
    double g1 = p[0]  * TMath::Gaus(xx, p[1],  p[2],  true);
    double g2 = p[3]  * TMath::Gaus(xx, p[4],  p[5],  true);
    double g3 = p[6]  * TMath::Gaus(xx, p[7],  p[8],  true);

    // --- spline ---
    double splineVal = 0.0;
    if (gSplinePS) splineVal = p[9] * gSplinePS->Eval(xx);

    // --- nuevos picos ---
    double g4 = p[10] * TMath::Gaus(xx, p[11], p[12], true); // 3.2 MeV
    double g5 = p[13] * TMath::Gaus(xx, p[14], p[15], true); // 5.4 MeV

    return g1 + g2 + g3 + splineVal + g4 + g5;
}

void fit_7Li_dp() {
    // Get file and Ex histo for 8Li excitation energy spectrum
    TFile *file = TFile::Open("../../PostAnalysis/Outputs/tree_ex_7Li_d_p.root");
    TH1D *hEx_7Li_dp = (TH1D *)file->Get("hExSil");

    // Get file and Ex for the PS simulated
    TFile *filePS = TFile::Open("./Inputs/1nPS_7Li.root");
    TH1D *hPS_1n = (TH1D *)filePS->Get("hEx");

    // spline global
    gSplinePS = new TSpline3(hPS_1n);

    // TF1 con 16 parámetros (3*3 + 1 + 2*3)
    TF1 *fTotal = new TF1("fTotal", fitFunc, -2, 6, 16);

    // nombres de parámetros: hay que asignarlos de a uno
    fTotal->SetParName(0,"A1");
    fTotal->SetParName(1,"mean1");
    fTotal->SetParName(2,"sigma1");
    fTotal->SetParName(3,"A2");
    fTotal->SetParName(4,"mean2");
    fTotal->SetParName(5,"sigma2");
    fTotal->SetParName(6,"A3");
    fTotal->SetParName(7,"mean3");
    fTotal->SetParName(8,"sigma3");
    fTotal->SetParName(9,"AmpSpline");
    fTotal->SetParName(10,"A4");
    fTotal->SetParName(11,"mean4");
    fTotal->SetParName(12,"sigma4");
    fTotal->SetParName(13,"A5");
    fTotal->SetParName(14,"mean5");
    fTotal->SetParName(15,"sigma5");

    // parámetros iniciales: usar un array
    double initPars[16] = {
        700, 0,   0.08,  // g1
        30, 1,   0.08,  // g2
        60, 2,   0.08,  // g3
        100,              // spline
        20, 3.2, 0.08,  // g4
        15, 5.4, 0.08   // g5
    };
    fTotal->SetParameters(initPars);

    // Set limits for parameters
    for (int i = 0; i < 16; i++) 
    {
        if(i == 1)
            fTotal->SetParLimits(i, -0.5, 0.5);
        else if(i == 4)
            fTotal->SetParLimits(i, 0.5, 1.5);
        else if(i == 7)
            fTotal->SetParLimits(i, 1.5, 2.5);
        else if(i == 11)
            fTotal->SetParLimits(i, 2.5, 3.9);
        else if(i == 14)
            fTotal->SetParLimits(i, 4.5, 6.0);
        else
            fTotal->SetParLimits(i, 0.01, 1e6); // mínimo=0, máximo grande
    }

    // fit
    hEx_7Li_dp->Fit(fTotal,"R");

    // Plot resultado
    TCanvas *c1 = new TCanvas("c1", "Ex 7Li(d,p)", 900, 700);
    hEx_7Li_dp->Draw();
    fTotal->Draw("same");

    // spline check
    TCanvas *c2 = new TCanvas("c2", "PS 1n 7Li", 800, 600);
    hPS_1n->Draw();
    gSplinePS->Draw("same");
}
