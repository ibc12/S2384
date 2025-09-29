#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TSpline.h"
#include "TMath.h"
#include "TLegend.h"
#include "TStyle.h"
#include <vector>

TSpline3 *gSplinePS = nullptr;
int gNumGauss = 0;
int gNumVoigt = 0;
int gNumPS = 0;

// --- Función genérica para calcular fit ---
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

    // Phase-space (spline)
    for(int i=0;i<gNumPS;i++){
        int idx = gNumGauss*3 + gNumVoigt*4 + i;
        if(gSplinePS) val += p[idx]*gSplinePS->Eval(xx);
    }

    return val;
}

// --- Función para inicializar parámetros ---
void InitFitParameters(TF1* fTotal,
                       const std::vector<double>& gaussInitMean,
                       const std::vector<double>& gaussInitSigma,
                       const std::vector<double>& voigtInitMean,
                       const std::vector<double>& voigtInitSigma,
                       const std::vector<double>& voigtInitGamma,
                       const std::vector<double>& ampInitGauss,
                       const std::vector<double>& ampInitVoigt,
                       const std::vector<double>& ampInitPS,
                       const std::vector<double>& gaussMeanMin,
                       const std::vector<double>& gaussMeanMax,
                       const std::vector<double>& gaussSigmaMin,
                       const std::vector<double>& gaussSigmaMax,
                       const std::vector<double>& voigtMeanMin,
                       const std::vector<double>& voigtMeanMax,
                       const std::vector<double>& voigtSigmaMin,
                       const std::vector<double>& voigtSigmaMax,
                       const std::vector<double>& voigtGammaMin,
                       const std::vector<double>& voigtGammaMax,
                       const std::vector<double>& ampGaussMin,
                       const std::vector<double>& ampGaussMax,
                       const std::vector<double>& ampVoigtMin,
                       const std::vector<double>& ampVoigtMax,
                       const std::vector<double>& ampPSMin,
                       const std::vector<double>& ampPSMax)
{
    int idx = 0;

    // Gaussianas
    for(size_t i=0;i<gaussInitMean.size();i++){
        fTotal->SetParName(idx,Form("A_G%d",int(i+1)));
        fTotal->SetParameter(idx,ampInitGauss[i]);
        fTotal->SetParLimits(idx,ampGaussMin[i],ampGaussMax[i]);
        idx++;

        fTotal->SetParName(idx,Form("mean_G%d",int(i+1)));
        fTotal->SetParameter(idx,gaussInitMean[i]);
        fTotal->SetParLimits(idx,gaussMeanMin[i],gaussMeanMax[i]);
        idx++;

        fTotal->SetParName(idx,Form("sigma_G%d",int(i+1)));
        fTotal->SetParameter(idx,gaussInitSigma[i]);
        fTotal->SetParLimits(idx,gaussSigmaMin[i],gaussSigmaMax[i]);
        idx++;
    }

    // Voigts
    for(size_t i=0;i<voigtInitMean.size();i++){
        fTotal->SetParName(idx,Form("A_V%d",int(i+1)));
        fTotal->SetParameter(idx,ampInitVoigt[i]);
        fTotal->SetParLimits(idx,ampVoigtMin[i],ampVoigtMax[i]);
        idx++;

        fTotal->SetParName(idx,Form("mean_V%d",int(i+1)));
        fTotal->SetParameter(idx,voigtInitMean[i]);
        fTotal->SetParLimits(idx,voigtMeanMin[i],voigtMeanMax[i]);
        idx++;

        fTotal->SetParName(idx,Form("sigma_V%d",int(i+1)));
        fTotal->SetParameter(idx,voigtInitSigma[i]);
        fTotal->SetParLimits(idx,voigtSigmaMin[i],voigtSigmaMax[i]);
        idx++;

        fTotal->SetParName(idx,Form("gamma_V%d",int(i+1)));
        fTotal->SetParameter(idx,voigtInitGamma[i]);
        fTotal->SetParLimits(idx,voigtGammaMin[i],voigtGammaMax[i]);
        idx++;
    }

    // Phase-space
    for(size_t i=0;i<ampInitPS.size();i++){
        fTotal->SetParName(idx,Form("PS%d",int(i+1)));
        fTotal->SetParameter(idx,ampInitPS[i]);
        fTotal->SetParLimits(idx,ampPSMin[i],ampPSMax[i]);
        idx++;
    }
}

void DrawFitComponentsGeneric(TH1D *histo, TF1 *fTotal, double fLim) {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    TCanvas *cComp = new TCanvas("cComp","Componentes del ajuste",900,700);
    histo->Draw();

    TLegend *leg = new TLegend(0.60,0.60,0.88,0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    Color_t colors[] = {kAzure+1, kOrange+7, kViolet+6, kTeal+3, kPink+9, kGreen+2, kRed+2};

    int paramIdx = 0;

    // Gaussianas
    for(int i=0;i<gNumGauss;i++){
        TF1 *g = new TF1(Form("g%d",i+1),
                         "[0]*TMath::Gaus(x,[1],[2],true)", -2, fLim);
        double A    = fTotal->GetParameter(paramIdx);
        double mean = fTotal->GetParameter(paramIdx+1);
        double sigma= fTotal->GetParameter(paramIdx+2);
        g->SetParameters(A, mean, sigma);
        g->SetLineColor(colors[i % 7]);
        g->Draw("same");
        leg->AddEntry(g, Form("g%d (%.2f MeV, #sigma=%.0f keV)", i+1, mean, sigma*1000.0), "l");
        paramIdx += 3;
    }

    // Voigts
    for(int i=0;i<gNumVoigt;i++){
        int idx = paramIdx;
        double A     = fTotal->GetParameter(idx);
        double mean  = fTotal->GetParameter(idx+1);
        double sigma = fTotal->GetParameter(idx+2);
        double gamma = fTotal->GetParameter(idx+3);

        // TF1 directo, sin lambda
        TF1 *v = new TF1(Form("v%d",i+1),
                         "[0]*TMath::Voigt(x-[1],[2],[3])", -2, fLim);
        v->SetParameters(A, mean, sigma, gamma);
        v->SetLineColor(colors[(i+gNumGauss) % 7]);
        v->Draw("same");
        leg->AddEntry(v, Form("v%d (%.2f MeV, #sigma=%.0f keV, #Gamma=%.0f keV)", i+1, mean, sigma*1000.0, gamma*1000.0), "l");

        paramIdx += 4;
    }

    // Phase-space (splines)
    for(int i=0;i<gNumPS;i++){
        double amp = fTotal->GetParameter(paramIdx);
        TF1 *ps = new TF1(Form("PS%d",i+1),
                          [amp](double *x,double *p){ return amp*gSplinePS->Eval(x[0]); }, -2,fLim,0);
        ps->SetLineColor(kGray+2);
        ps->SetLineStyle(2);
        ps->Draw("same");
        leg->AddEntry(ps, Form("PS%d", i+1), "l");
        paramIdx += 1;
    }

    leg->Draw();
}


// --- Macro principal ---
void fit_dd()
{
    // --- Configuración ---
    const char* histoFile = "../../PostAnalysis/Outputs/tree_ex_7Li_d_d.root";
    const char* histoName = "hExSil";
    const char* psFile    {};
    double fLim = 8.5;
    int rebinFactor = 2;

    // --- Números de picos ---
    gNumGauss = 1;
    gNumVoigt = 0;
    gNumPS    = 0;

    // --- Abrir histograma ---
    TFile *file = TFile::Open(histoFile);
    TH1D *h = (TH1D*)file->Get(histoName);
    if(rebinFactor>1) 
    {
        h->Rebin(rebinFactor);
        TString title {Form("Ex with silicons;E_{x} [MeV];Counts / (%i keV)", 75*rebinFactor)};
        h->SetTitle(title);
    }
    // --- Spline ---
    if(psFile){
        TFile *fPS = TFile::Open(psFile);
        TH1D *hPS = (TH1D*)fPS->Get("hEx");
        gSplinePS = new TSpline3(hPS);
    }

    // --- TF1 total ---
    int nParTotal = gNumGauss*3 + gNumVoigt*4 + gNumPS;
    TF1 *fTotal = new TF1("fTotal",fitFuncGeneric,-2,fLim,nParTotal);

    // --- Declaración de vectores de parámetros iniciales ---
    std::vector<double> gaussMean   = {0};
    std::vector<double> gaussSigma  = {0.08};
    std::vector<double> voigtMean   = {};
    std::vector<double> voigtSigma  = {};
    std::vector<double> voigtGamma  = {};
    std::vector<double> ampGauss    = {700};
    std::vector<double> ampVoigt    = {};
    std::vector<double> ampPS       = {};

    // --- Declaración de vectores de límites ---
    std::vector<double> gaussMeanMin = {-0.5};
    std::vector<double> gaussMeanMax = {0.5};
    std::vector<double> gaussSigmaMin= {0.05};
    std::vector<double> gaussSigmaMax= {0.5};

    std::vector<double> voigtMeanMin = {};
    std::vector<double> voigtMeanMax = {};
    std::vector<double> voigtSigmaMin= {};
    std::vector<double> voigtSigmaMax= {};
    std::vector<double> voigtGammaMin= {};
    std::vector<double> voigtGammaMax= {};

    std::vector<double> ampGaussMin = {0.1};
    std::vector<double> ampGaussMax = {10000};
    std::vector<double> ampVoigtMin = {};
    std::vector<double> ampVoigtMax = {};
    std::vector<double> ampPSMin    = {};
    std::vector<double> ampPSMax    = {};

    // --- Inicializar parámetros y límites ---
    InitFitParameters(fTotal,
                      gaussMean, gaussSigma,
                      voigtMean, voigtSigma, voigtGamma,
                      ampGauss, ampVoigt, ampPS,
                      gaussMeanMin, gaussMeanMax,
                      gaussSigmaMin, gaussSigmaMax,
                      voigtMeanMin, voigtMeanMax,
                      voigtSigmaMin, voigtSigmaMax,
                      voigtGammaMin, voigtGammaMax,
                      ampGaussMin, ampGaussMax,
                      ampVoigtMin, ampVoigtMax,
                      ampPSMin, ampPSMax);

    // --- Fit ---
    h->Fit(fTotal,"R");

    // --- Canvas total ---
    TCanvas *c1 = new TCanvas("c1","Fit total",900,700);
    h->Draw();
    fTotal->Draw("same");

    // --- Canvas componentes ---
    // (Puedes usar la función DrawFitComponentsGeneric si la tienes definida)
    DrawFitComponentsGeneric(h, fTotal, fLim);

    // Guardar el histograma y la función
    TFile *fout = new TFile("./Outputs/ExFitTotal.root","RECREATE");
    h->Write("hExTotal");
    fTotal->Write("fExTotal");
    if(gSplinePS)
        gSplinePS->Write("splinePS");
    fout->Close();
}
