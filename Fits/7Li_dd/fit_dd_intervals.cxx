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

#include "ROOT/RDataFrame.hxx"

#include "ActMergerData.h"

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

// fAmps debe tener tamaño: gNumGauss + gNumVoigt + gNumPS
void FixParametersFromReference(TF1* fTarget, TF1* fRef,
                                const std::vector<double>& fAmps = {}) {
    int paramIdx = 0;

    // --- Gaussianas ---
    for(int i = 0; i < gNumGauss; i++) {
        // Nombre y parámetros de referencia
        fTarget->SetParName(paramIdx, fRef->GetParName(paramIdx));       // amplitud
        fTarget->SetParName(paramIdx+1, fRef->GetParName(paramIdx+1));   // posición
        fTarget->SetParName(paramIdx+2, fRef->GetParName(paramIdx+2));   // sigma

        // Posición y sigma fijadas
        fTarget->SetParameter(paramIdx+1, fRef->GetParameter(paramIdx+1));
        fTarget->FixParameter(paramIdx+1, fRef->GetParameter(paramIdx+1));

        fTarget->SetParameter(paramIdx+2, fRef->GetParameter(paramIdx+2));
        fTarget->FixParameter(paramIdx+2, fRef->GetParameter(paramIdx+2));

        // Amplitud inicial
        if(!fAmps.empty() && i < fAmps.size()) fTarget->SetParameter(paramIdx, fAmps[i]);

        // Liberar solo amplitud
        fTarget->ReleaseParameter(paramIdx);

        paramIdx += 3;
    }

    // --- Voigts ---
    for(int i = 0; i < gNumVoigt; i++) {
        fTarget->SetParName(paramIdx, fRef->GetParName(paramIdx));       // amplitud
        for(int j = 1; j < 4; j++) fTarget->SetParName(paramIdx+j, fRef->GetParName(paramIdx+j));

        // Fijar posición y anchos
        for(int j = 1; j < 4; j++) {
            fTarget->SetParameter(paramIdx+j, fRef->GetParameter(paramIdx+j));
            fTarget->FixParameter(paramIdx+j, fRef->GetParameter(paramIdx+j));
        }

        // Amplitud inicial
        if(!fAmps.empty() && (gNumGauss + i) < fAmps.size()) 
            fTarget->SetParameter(paramIdx, fAmps[gNumGauss+i]);

        // Liberar solo amplitud
        fTarget->ReleaseParameter(paramIdx);

        paramIdx += 4;
    }

    // --- Phase-space ---
    for(int i = 0; i < gNumPS; i++) {
        if(!fAmps.empty() && (gNumGauss+gNumVoigt+i) < fAmps.size()) 
            fTarget->SetParameter(paramIdx, fAmps[gNumGauss+gNumVoigt+i]);

        fTarget->ReleaseParameter(paramIdx); // liberar amplitud
        paramIdx += 1;
    }
}



void GetIntervalHists(TFile *fData, std::vector<TH1D*> &hists,
                      TH1D* hRef,
                      const std::vector<double>& minAngle,
                      const std::vector<double>& maxAngle) {
    // Crear dataframe
    ROOT::RDataFrame df("Final_Tree", fData);
    auto def = df.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled(); }, {"MergerData"});

    if (!hRef) {
        std::cout << "Error: histograma de referencia nulo." << std::endl;
        return;
    }

    // Usar el binning exacto de hRef
    int nbins   = hRef->GetNbinsX();
    double xmin = hRef->GetXaxis()->GetXmin();
    double xmax = hRef->GetXaxis()->GetXmax();

    // Loop sobre intervalos angulares
    for (size_t i = 0; i < minAngle.size(); ++i) {
        auto df_filt = def.Filter(Form("ThetaCM >= %f && ThetaCM < %f",
                                       minAngle[i], maxAngle[i]));

        std::string hname = Form("hEx_theta_%zu", i);
        auto h_res = df_filt.Histo1D(
            {hname.c_str(),
             Form("Ex %.0f-%.0f deg", minAngle[i], maxAngle[i]),
             nbins, xmin, xmax},  // <-- mismo binning que hRef
            "Ex");

        // Obtener el TH1D* y clonarlo para que sobreviva
        TH1D* hclone = (TH1D*)h_res.GetPtr()->Clone(hname.c_str());
        hclone->SetDirectory(0);
        hists.push_back(hclone);
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
void fit_dd_intervals() {
    // --- Configuración ---
    gNumGauss = 1;  // deben coincidir con la macro global
    gNumVoigt = 0;
    gNumPS    = 0;
    double fLim = 8.5;

    // === 1. Cargar ajuste global ===
    TFile *fRef = TFile::Open("./Outputs/ExFitTotal.root");
    TH1D *hRef = (TH1D*)fRef->Get("hExTotal");
    TF1  *fRefFunc = (TF1*)fRef->Get("fExTotal");
    TSpline3 *splineRef = (TSpline3*)fRef->Get("splinePS");
    gSplinePS = splineRef; // asignar spline global

    // === 2. Cargar histos de intervalos ===
    // Angular intervals
    std::vector<double> minAngle {10, 25, 40, 55, 70};
    std::vector<double> maxAngle {25, 40, 55, 70, 85};
    // File and getter
    TFile *fData = TFile::Open("../../PostAnalysis/Outputs/tree_ex_7Li_d_d.root"); 
    std::vector<TH1D*> hists;
    GetIntervalHists(fData, hists, hRef, minAngle, maxAngle);

    TH1D *h = (TH1D*)fData->Get("hExSil");

    TFile *fParams = TFile::Open("./Outputs/ExFitTotal.root");
    TF1  *fExTotal = (TF1*)fParams->Get("fExTotal");

    if(!h || !fExTotal) {
        std::cerr << "Error: No se pudo cargar el histograma o la función de referencia." << std::endl;
        return;
    }

    std::vector<double> AmpIni = {20, 10}; // Initial amplitudes for fit
    // === 3. Loop ===
    for (size_t i=0;i<hists.size();i++) {
        TH1D *h = hists[i];

        int nParTotal = gNumGauss*3 + gNumVoigt*4 + gNumPS;
        TF1 *fAmp = new TF1(Form("fAmp_theta%zu",i),
                            fitFuncGeneric,-2,fLim,nParTotal);

        FixParametersFromReference(fAmp,fRefFunc, AmpIni);

        h->Fit(fAmp,"R");

        // Dibujo
        DrawFitComponentsGeneric(h,fAmp,fLim);
    }
}
