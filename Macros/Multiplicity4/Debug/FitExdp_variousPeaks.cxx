#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActParticle.h"
#include "ActSRIM.h"
#include "ActSilData.h"
#include "ActSilMatrix.h"
#include "ActSilSpecs.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TF1.h"
#include "TH2.h"
#include "TH2D.h"
#include "TPaveText.h"
#include "TString.h"

#include <fstream>
#include <map>
#include <string>

#include "../../../Fits//Histos.h"
#include "../../../PostAnalysis/HistConfig.h"
#include "../../../PrettyStyle.C"
#include "../Utils.h"


void FitExdp_variousPeaks()
{
    PrettyStyle(false, false);
    // ROOT::EnableImplicitMT();

    // Analysis
    ROOT::RDataFrame df {"Final_Tree", "../Outputs/ExM4_7Li_d_p.root"};
    auto def = df.Filter([](bool StoppedInActar) { return StoppedInActar == false; },
                         {"StoppedInACTAR"}); // only silicons, == false is for L1 events

    // Histogram of Ex
    auto hEx {def.Histo1D(S2384Fit::Exdp_7Li, "Ex")};
    hEx->SetAxisRange(4, 10);
    hEx->SetStats(0);

    auto hExRebined2 {def.Histo1D(S2384Fit::Exdp_7Li, "Ex")};
    hExRebined2->Rebin(2);
    hExRebined2->SetAxisRange(4, 10);
    hExRebined2->GetYaxis()->SetTitle("Counts / 300 keV");
    hExRebined2->SetStats(0);

    auto hExRebined3 {def.Histo1D(S2384Fit::Exdp_7Li, "Ex")};
    hExRebined3->Rebin(3);
    hExRebined3->SetAxisRange(4, 10);
    hExRebined3->GetYaxis()->SetTitle("Counts / 450 keV");
    hExRebined3->SetStats(0);

    auto seedAmplitude = [](TH1* h, double meanGuess, double gammaGuess)
    {
        double H = h->GetBinContent(h->FindBin(meanGuess));
        return H * TMath::Pi() * gammaGuess / 2.;
    };

    auto makeSingleBW = [](const char* name)
    {
        auto* f = new TF1(name, "[0]*TMath::BreitWigner(x,[1],[2])", 5.5, 9.5);
        f->SetParNames("A", "E_{r}", "#Gamma");
        f->SetParLimits(2, 0.01, 1.0);
        f->SetLineColor(kRed);
        return f;
    };

    auto makeDoubleBW = [](const char* name)
    {
        auto* f = new TF1(name, "[0]*TMath::BreitWigner(x,[1],[2]) + [3]*TMath::BreitWigner(x,[4],[5])", 5.5, 9.5);
        f->SetParNames("A_{1}", "E_{1}", "#Gamma_{1}", "A_{2}", "E_{2}", "#Gamma_{2}");
        f->SetParLimits(1, 6.0, 7.0);
        f->SetParLimits(2, 0.01, 1.0);
        f->SetParLimits(4, 7.0, 8.0);
        f->SetParLimits(5, 0.01, 1.0);
        f->SetLineColor(kRed);
        return f;
    };

    auto* fS1 {makeSingleBW("fBW_hEx")};
    auto* fS2 {makeSingleBW("fBW_hExRebined2")};
    auto* fS3 {makeSingleBW("fBW_hExRebined3")};

    auto* f2_1 {makeDoubleBW("f2BW_hEx")};
    auto* f2_2 {makeDoubleBW("f2BW_hExRebined2")};
    auto* f2_3 {makeDoubleBW("f2BW_hExRebined3")};

    fS1->SetParameters(seedAmplitude(hEx.GetPtr(), 7.0, 0.5), 7.0, 0.5);
    fS2->SetParameters(seedAmplitude(hExRebined2.GetPtr(), 7.0, 0.5), 7.0, 0.5);
    fS3->SetParameters(seedAmplitude(hExRebined3.GetPtr(), 7.0, 0.5), 7.0, 0.5);

    f2_1->SetParameters(seedAmplitude(hEx.GetPtr(), 6.5, 0.3), 6.5, 0.3, seedAmplitude(hEx.GetPtr(), 7.5, 0.3), 7.5,
                        0.3);
    f2_2->SetParameters(seedAmplitude(hExRebined2.GetPtr(), 6.5, 0.3), 6.5, 0.3,
                        seedAmplitude(hExRebined2.GetPtr(), 7.5, 0.3), 7.5, 0.3);
    f2_3->SetParameters(seedAmplitude(hExRebined3.GetPtr(), 6.5, 0.3), 6.5, 0.3,
                        seedAmplitude(hExRebined3.GetPtr(), 7.5, 0.3), 7.5, 0.3);

    // SIEMPRE en la mitad izquierda del pad: x1=0.13, x2=0.48
    auto drawFitText = [](TF1* f, double y1, double y2)
    {
        auto* pt = new TPaveText(0.13, y1, 0.48, y2, "NDC");
        pt->SetFillColor(0);
        pt->SetFillStyle(0);
        pt->SetBorderSize(0);
        pt->SetTextAlign(12);
        pt->SetTextFont(42);
        pt->SetTextSize(0.032);

        double chi2 {f->GetChisquare()};
        int ndf {f->GetNDF()};
        pt->AddText(Form("#chi^{2}/ndf = %.2f / %d = %.2f", chi2, ndf, ndf > 0 ? chi2 / ndf : 0.));
        for(int i = 0; i < f->GetNpar(); i++)
        {
            pt->AddText(Form("%s = %.3f #pm %.3f", f->GetParName(i), f->GetParameter(i), f->GetParError(i)));
        }
        pt->Draw();
        return pt;
    };

    auto* c = new TCanvas("cExdp_variousPeaks", "Excitation energy Multiplicity 4 - L1", 1200, 800);
    c->DivideSquare(6);

    // Fila 1: fit simple (4 líneas: chi2 + 3 params) -> caja más baja
    c->cd(1);
    hEx->Fit(fS1, "R");
    hEx->DrawClone();
    drawFitText(fS1, 0.65, 0.89);

    c->cd(2);
    hExRebined2->Fit(fS2, "R");
    hExRebined2->DrawClone();
    drawFitText(fS2, 0.65, 0.89);

    c->cd(3);
    hExRebined3->Fit(fS3, "R");
    hExRebined3->DrawClone();
    drawFitText(fS3, 0.65, 0.89);

    // Fila 2: fit doble (7 líneas: chi2 + 6 params) -> caja más alta, MISMO x1/x2 (izquierda)
    c->cd(4);
    hEx->Fit(f2_1, "R");
    hEx->DrawClone();
    drawFitText(f2_1, 0.40, 0.89);

    c->cd(5);
    hExRebined2->Fit(f2_2, "R");
    hExRebined2->DrawClone();
    drawFitText(f2_2, 0.40, 0.89);

    c->cd(6);
    hExRebined3->Fit(f2_3, "R");
    hExRebined3->DrawClone();
    drawFitText(f2_3, 0.40, 0.89);

    c->cd();
    c->Update();
}