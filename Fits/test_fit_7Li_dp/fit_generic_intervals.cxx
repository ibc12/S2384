#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TSpline.h"
#include "TMath.h"
#include "TLegend.h"
#include "TStyle.h"
#include <vector>
#include <string>

TSpline3 *gSplinePS = nullptr;
int gNumGauss = 0;
int gNumVoigt = 0;
int gNumPS    = 0;

// --- Función genérica de fit ---
double fitFuncGeneric(double *x, double *p) {
    double xx = x[0];
    double val = 0;

    // Gaussianas
    for(int i=0;i<gNumGauss;i++){
        int idx = i*3;
        val += p[idx]*TMath::Gaus(xx,p[idx+1],p[idx+2],true);
    }

    // Voigts
    for(int i=0;i<gNumVoigt;i++){
        int idx = gNumGauss*3 + i*4;
        val += p[idx]*TMath::Voigt(xx-p[idx+1],p[idx+2],p[idx+3]);
    }

    // Phase-space
    for(int i=0;i<gNumPS;i++){
        int idx = gNumGauss*3 + gNumVoigt*4 + i;
        if(gSplinePS) val += p[idx]*gSplinePS->Eval(xx);
    }

    return val;
}

// --- Fijar todos los parámetros al valor de referencia ---
void FixParametersFromReference(TF1* fTarget, TF1* fRef) {
    for (int i=0;i<fRef->GetNpar();i++) {
        fTarget->SetParName(i, fRef->GetParName(i));
        fTarget->SetParameter(i, fRef->GetParameter(i));
        fTarget->FixParameter(i, fRef->GetParameter(i));
    }

    // liberar amplitudes
    int idx = 0;
    for (int g=0; g<gNumGauss; g++) {
        fTarget->ReleaseParameter(idx);
        idx += 3;
    }
    for (int v=0; v<gNumVoigt; v++) {
        fTarget->ReleaseParameter(idx);
        idx += 4;
    }
    for (int ps=0; ps<gNumPS; ps++) {
        fTarget->ReleaseParameter(idx);
        idx += 1;
    }
}

// --- Dibujar componentes ---
void DrawFitComponentsGeneric(TH1D *histo, TF1 *fTotal, double fLim) {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    TCanvas *cComp = new TCanvas(Form("cComp_%s",histo->GetName()),
                                 "Componentes del ajuste",900,700);
    histo->Draw();

    TLegend *leg = new TLegend(0.60,0.60,0.88,0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    Color_t colors[] = {kAzure+1, kOrange+7, kViolet+6, kTeal+3,
                        kPink+9, kGreen+2, kRed+2};

    int paramIdx = 0;

    // Gaussianas
    for(int i=0;i<gNumGauss;i++){
        TF1 *g = new TF1(Form("g%d_%s",i+1,histo->GetName()),
                         "[0]*TMath::Gaus(x,[1],[2],true)",
                         -2, fLim);
        g->SetParameters(fTotal->GetParameter(paramIdx),
                         fTotal->GetParameter(paramIdx+1),
                         fTotal->GetParameter(paramIdx+2));
        g->SetLineColor(colors[i % 7]);
        g->Draw("same");
        leg->AddEntry(g, Form("g%d (%.2f MeV, #sigma=%.0f keV)",
                              i+1,
                              fTotal->GetParameter(paramIdx+1),
                              fTotal->GetParameter(paramIdx+2)*1000.0), "l");
        paramIdx += 3;
    }

    // Voigts
    for(int i=0;i<gNumVoigt;i++){
        TF1 *v = new TF1(Form("v%d_%s",i+1,histo->GetName()),
                         "[0]*TMath::Voigt(x-[1],[2],[3])",
                         -2, fLim);
        v->SetParameters(fTotal->GetParameter(paramIdx),
                         fTotal->GetParameter(paramIdx+1),
                         fTotal->GetParameter(paramIdx+2),
                         fTotal->GetParameter(paramIdx+3));
        v->SetLineColor(colors[(i+gNumGauss) % 7]);
        v->Draw("same");
        leg->AddEntry(v, Form("v%d (%.2f MeV, #sigma=%.0f keV, #Gamma=%.0f keV)",
                              i+1,
                              fTotal->GetParameter(paramIdx+1),
                              fTotal->GetParameter(paramIdx+2)*1000.0,
                              fTotal->GetParameter(paramIdx+3)*1000.0), "l");
        paramIdx += 4;
    }

    // PS
    for(int i=0;i<gNumPS;i++){
        double amp = fTotal->GetParameter(paramIdx);
        TF1 *ps = new TF1(Form("PS%d_%s",i+1,histo->GetName()),
                          [amp](double *x,double*){ return amp*gSplinePS->Eval(x[0]); },
                          -2, fLim, 0);
        ps->SetLineColor(kGray+2);
        ps->SetLineStyle(2);
        ps->Draw("same");
        leg->AddEntry(ps, Form("PS%d", i+1), "l");
        paramIdx++;
    }

    leg->Draw();
}

// --- Macro principal para intervalos ---
void fit_generic_intervals() {
    // --- Configuración ---
    gNumGauss = 3;  // deben coincidir con la macro global
    gNumVoigt = 4;
    gNumPS    = 1;
    double fLim = 8.5;

    // === 1. Cargar ajuste global ===
    TFile *fRef = TFile::Open("./Outputs/ExFitTotal.root");
    TH1D *hRef = (TH1D*)fRef->Get("hExTotal");
    TF1  *fRefFunc = (TF1*)fRef->Get("fExTotal");

    // === 2. Cargar histos de intervalos ===
    TFile *fIn = TFile::Open("../../PostAnalysis/Outputs/tree_ex_7Li_d_p.root"); 
    std::vector<TH1D*> hists;
    for (int i=0;i<1;i++) {
        TH1D *h = (TH1D*)fIn->Get("hExSil"); 
        if(h) hists.push_back(h);
    }

    TFile *fData = TFile::Open("../../PostAnalysis/Outputs/tree_ex_7Li_d_p.root"); 
    TH1D *h = (TH1D*)fData->Get("hExSil");

    TFile *fParams = TFile::Open("./Outputs/ExFitTotal.root");
    TF1  *fExTotal = (TF1*)fParams->Get("fExTotal");

    if(!h || !fExTotal) {
        std::cerr << "Error: No se pudo cargar el histograma o la función de referencia." << std::endl;
        return;
    }

    // === 3. Loop ===
    for (size_t i=0;i<hists.size();i++) {
        TH1D *h = hists[i];

        int nParTotal = gNumGauss*3 + gNumVoigt*4 + gNumPS;
        TF1 *fAmp = new TF1(Form("fAmp_theta%zu",i),
                            fitFuncGeneric,-2,fLim,nParTotal);

        FixParametersFromReference(fAmp,fRefFunc);

        h->Fit(fAmp,"R");

        // Dibujo
        DrawFitComponentsGeneric(h,fAmp,fLim);
    }
}
