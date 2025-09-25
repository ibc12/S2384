#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TSpline.h"
#include "TMath.h"
#include "TLegend.h"
#include "TStyle.h"

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

void DrawFitComponents(TH1D *histo, TF1 *fTotal, double fLim, TSpline3 *spline = nullptr) {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0); // quitar también la caja del fit
    
    TCanvas *cComp = new TCanvas("cComp","Componentes del ajuste",900,700);
    histo->Draw();

    // Colores para gaussianas
    Color_t colors[] = {kAzure+1, kOrange+7, kViolet+6, kTeal+3, kPink+9};

    // Crear la leyenda en la parte superior derecha
    TLegend *leg = new TLegend(0.60,0.60,0.88,0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    // Dibujar las 5 gaussianas
    const char* labels[5] = {"g.s.", "1st", "2nd", "3rd", "4th"};
    for (int j = 0; j < 5; j++) {
        int idx = j*3;
        if (j >= 3) idx += 1; // saltar AmpSpline

        double mean  = fTotal->GetParameter(idx+1); // en MeV
        double sigma = fTotal->GetParameter(idx+2); // en MeV

        TF1 *g = new TF1(Form("g%d",j+1),
                         "[0]*TMath::Gaus(x,[1],[2],true)", -2, fLim);
        g->SetParameters(fTotal->GetParameter(idx),
                         mean, sigma);
        g->SetLineColor(colors[j]);
        g->Draw("same");

        // Texto con nombre + mean en MeV + sigma en keV
        TString entry = Form("%s (%.2f MeV, #sigma=%.0f keV)", 
                             labels[j], mean, sigma*1000.0);
        leg->AddEntry(g, entry, "l");
    }

    // spline (si no está fijada a 0)
    if (spline && fTotal->GetParameter(9) != 0.0) {
        TF1 *fspline = new TF1("fspline",
            [spline](double *x,double *p){ return p[0]*spline->Eval(x[0]); },
            -2,fLim,1);
        fspline->SetParameter(0, fTotal->GetParameter(9));
        fspline->SetLineColor(kGray+2);
        fspline->SetLineStyle(2);
        fspline->Draw("same");

        leg->AddEntry(fspline, "PS (S_{n} = 2.04 MeV)", "l");
    }

    leg->Draw();
}

void fit_7Li_dp() {
    // Get file and Ex histo
    TFile *file = TFile::Open("../../PostAnalysis/Outputs/tree_ex_7Li_d_p.root");
    TH1D *hEx_7Li_dp = (TH1D *)file->Get("hExSil");
    bool rebin {true};
    if (rebin)
    {
        int rebinFactor = 2;
        hEx_7Li_dp->Rebin(2);
        TString title {Form("Ex with silicons;E_{x} [MeV];Counts / (%i keV)", 75*rebinFactor)};
        hEx_7Li_dp->SetTitle(title);
    }
    

    // Get PS file
    TFile *filePS = TFile::Open("./Inputs/1nPS_7Li_latSil.root");
    TH1D *hPS_1n = (TH1D *)filePS->Get("hEx");

    // spline global
    gSplinePS = new TSpline3(hPS_1n);

    // TF1 con 16 parámetros
    double fLim {6.5};
    TF1 *fTotal = new TF1("fTotal", fitFunc, -2, fLim, 16);

    // nombres de parámetros
    const char* names[16] = {
        "A1","mean1","sigma1","A2","mean2","sigma2","A3","mean3","sigma3",
        "AmpSpline","A4","mean4","sigma4","A5","mean5","sigma5"
    };
    for (int i=0;i<16;i++) fTotal->SetParName(i,names[i]);

    // parámetros iniciales
    double initPars[16] = {
        700, 0,   0.08,  // g1
        30, 1,   0.08,  // g2
        60, 2,   0.08,  // g3
        100,              // spline
        20, 3.2, 0.08,  // g4
        15, 5.4, 0.08   // g5
    };
    fTotal->SetParameters(initPars);

    // límites
    for (int i = 0; i < 16; i++) {
        //if (i == 9) continue;
        if(i == 1)       fTotal->SetParLimits(i, -0.5, 0.5); // g1 mean
        else if(i == 4)  fTotal->SetParLimits(i, 0.5, 1.5); // g2 mean
        else if(i == 7)  fTotal->SetParLimits(i, 2, 2.5); // g3 mean
        else if(i == 11) fTotal->SetParLimits(i, 2.7, 3.9); // g4 mean
        else if(i == 14) fTotal->SetParLimits(i, 4.5, 6.0); // g5 mean
        else if(i == 2 || i == 5 || i == 8 || i == 12 || i == 15) // sigmas
            fTotal->SetParLimits(i, 0.05, 0.5);
        else if(i == 13) fTotal->SetParLimits(i, 0.0, 30); // A5
        else
            fTotal->SetParLimits(i, 0.0, 100);
    }
    // fTotal->FixParameter(9, 0.);

    // fit
    hEx_7Li_dp->Fit(fTotal,"R");

    // Canvas con el ajuste total
    TCanvas *c1 = new TCanvas("c1", "Ex 7Li(d,p)", 900, 700);
    hEx_7Li_dp->Draw();
    fTotal->Draw("same");

    // Canvas con componentes
    DrawFitComponents(hEx_7Li_dp, fTotal, fLim, gSplinePS);

    // Canvas para revisar spline original
    TCanvas *c2 = new TCanvas("c2", "PS 1n 7Li", 800, 600);
    hPS_1n->Draw();
    gSplinePS->Draw("same");
}
