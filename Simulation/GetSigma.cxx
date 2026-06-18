#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

void GetSigma(const std::string& beam, const std::string& light, const std::string& target,
              const std::vector<double>& Exs, int neutronPS = 0, int protonPS = 0)
{
    const int N = static_cast<int>(Exs.size());
    if(N == 0)
    {
        std::cout << "[GetSigma] ERROR: Exs vector is empty.\n";
        return;
    }

    std::vector<double> vEx(N), vExErr(N, 0.0);
    std::vector<double> vSigma(N), vSigmaErr(N);

    gStyle->SetOptFit(1111);
    gStyle->SetOptStat("emr");

    std::string angStr {"_2-1AngStr"};

    std::vector<TH1D*> hVec;
    std::vector<TF1*> fVec;

    // ------------------------------------------------------------------ //
    //  Loop: open files, fill histos, fit
    // ------------------------------------------------------------------ //
    for(int i = 0; i < N; i++)
    {
        double ex = Exs[i];

        TString fname =
            TString::Format("./Outputs/%s/test_ang_straggling/%s_%s_TRIUMF_Eex_%.3f_nPS_%d_pPS_%d%s.root", beam.c_str(),
                            target.c_str(), light.c_str(), ex, neutronPS, protonPS, angStr.c_str());

        hVec.push_back(nullptr);
        fVec.push_back(nullptr);

        TFile* fin = TFile::Open(fname, "READ");
        if(!fin || fin->IsZombie())
        {
            std::cerr << "[GetSigma] WARNING: cannot open " << fname << " — skipping.\n";
            vEx[i] = ex;
            vSigma[i] = 0.0;
            vSigmaErr[i] = 0.0;
            continue;
        }

        TTree* tree = dynamic_cast<TTree*>(fin->Get("SimulationTTree"));
        if(!tree)
        {
            std::cerr << "[GetSigma] WARNING: TTree 'SimulationTTree' not found in " << fname << " — skipping.\n";
            fin->Close();
            vEx[i] = ex;
            vSigma[i] = 0.0;
            vSigmaErr[i] = 0.0;
            continue;
        }

        // ---------------------------------------------------------------- //
        //  Read Eex branch manually into a vector, then fill histogram
        //  (avoids the TTree::Draw + SetDirectory(nullptr) conflict)
        // ---------------------------------------------------------------- //
        double eexVal {};
        tree->SetBranchAddress("Eex", &eexVal);
        Long64_t nEntries = tree->GetEntries();

        double hMin = ex - 0.25;
        double hMax = ex + 0.25;

        TString hname = TString::Format("hEex_%d", i);
        TString htitle = TString::Format("E_{ex} = %.3f MeV;E_{ex} (MeV);Counts", ex);
        TH1D* hEex = new TH1D(hname, htitle, 200, hMin, hMax);
        hEex->SetDirectory(nullptr); // safe now because we fill manually

        for(Long64_t j = 0; j < nEntries; j++)
        {
            tree->GetEntry(j);
            hEex->Fill(eexVal);
        }

        if(hEex->GetEntries() == 0)
        {
            std::cerr << "[GetSigma] WARNING: histogram empty for Ex = " << ex << " — skipping.\n";
            fin->Close();
            delete hEex;
            vEx[i] = ex;
            vSigma[i] = 0.0;
            vSigmaErr[i] = 0.0;
            continue;
        }

        // ---------------------------------------------------------------- //
        //  Gaussian fit
        // ---------------------------------------------------------------- //
        TF1* fGaus = new TF1(TString::Format("fGaus_%d", i), "gaus", hMin, hMax);

        double initSigma = hEex->GetRMS() * 0.8;
        if(initSigma <= 0)
            initSigma = 0.05;

        fGaus->SetParameters(hEex->GetMaximum(), hEex->GetMean(), initSigma);
        fGaus->SetParLimits(2, 0.0, hMax - hMin);

        hEex->Fit(fGaus, "RQS");

        // Refit in ±2.5 sigma for better convergence
        double fMean = fGaus->GetParameter(1);
        double fSigma = std::abs(fGaus->GetParameter(2));
        if(fSigma > 0 && fSigma < (hMax - hMin))
        {
            fGaus->SetRange(std::max(fMean - 2.5 * fSigma, hMin), std::min(fMean + 2.5 * fSigma, hMax));
            hEex->Fit(fGaus, "RQS");
        }

        double sigma = std::abs(fGaus->GetParameter(2));
        double sigmaErr = fGaus->GetParError(2);

        std::cout << TString::Format("[GetSigma] Ex = %.3f MeV  ->  sigma = %.4f +/- %.4f MeV\n", ex, sigma, sigmaErr);

        vEx[i] = ex;
        vSigma[i] = sigma;
        vSigmaErr[i] = sigmaErr;

        hVec.back() = hEex;
        fVec.back() = fGaus;

        fin->Close();
        delete fin;
    }

    // ------------------------------------------------------------------ //
    //  Draw all histograms + fits on a properly divided canvas
    // ------------------------------------------------------------------ //
    int nCols = static_cast<int>(std::ceil(std::sqrt(static_cast<double>(N))));
    int nRows = static_cast<int>(std::ceil(static_cast<double>(N) / nCols));

    TCanvas* cFits = new TCanvas("cFits", "Eex Gaussian Fits", 400 * nCols, 350 * nRows);
    cFits->Divide(nCols, nRows);

    for(int i = 0; i < N; i++)
    {
        cFits->cd(i + 1);
        gPad->SetLeftMargin(0.13);

        if(!hVec[i])
            continue;

        hVec[i]->SetLineColor(kBlue + 1);
        hVec[i]->Draw("hist");

        if(fVec[i])
        {
            fVec[i]->SetLineColor(kRed);
            fVec[i]->SetLineWidth(2);
            fVec[i]->Draw("same");
        }
    }
    cFits->Update();

    // ------------------------------------------------------------------ //
    //  Build TGraphErrors and save
    // ------------------------------------------------------------------ //
    TGraphErrors* grSigma = new TGraphErrors(N, vEx.data(), vSigma.data(), vExErr.data(), vSigmaErr.data());
    grSigma->SetName("gsigmas");
    grSigma->SetTitle(
        TString::Format("Gaussian #sigma vs E_{ex}  [%s + %s #rightarrow %s + X];E_{ex} (MeV);#sigma (MeV)",
                        beam.c_str(), target.c_str(), light.c_str()));
    grSigma->SetMarkerStyle(20);
    grSigma->SetMarkerColor(kBlue + 1);
    grSigma->SetLineColor(kBlue + 1);

    TString outName = TString::Format("./Outputs/%s/test_ang_straggling/sigmas_%s_%s_%s%s.root", beam.c_str(),
                                      beam.c_str(), target.c_str(), light.c_str(), angStr.c_str());

    TFile* fout = TFile::Open(outName, "RECREATE");
    if(!fout || fout->IsZombie())
    {
        std::cerr << "[GetSigma] ERROR: cannot create output file " << outName << '\n';
        return;
    }
    grSigma->Write();
    cFits->Write();
    fout->Close();

    std::cout << "[GetSigma] Results saved to " << outName << '\n';

    TCanvas* cSigma = new TCanvas("cSigma", "Sigma vs Ex", 700, 500);
    grSigma->Draw("APE");
    cSigma->Update();
}