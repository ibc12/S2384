#include "ActDataManager.h"
#include "ActSilData.h"
#include "ActTPCData.h"
#include "ActModularData.h"
#include "ActModularData.h"
#include "ActMergerData.h"
#include "ActTPCData.h"

#include <ROOT/RDataFrame.hxx>

#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"
#include "TLegend.h"

void CFA_F2_comparison()
{
    // Get Sil data
    ActRoot::DataManager dataman{"../configs/data_all.conf", ActRoot::ModeType::EReadSilMod};
    auto chain{dataman.GetChain()};

    auto friend1{dataman.GetChain(ActRoot::ModeType::EMerge)};
    chain->AddFriend(friend1.get());

    // Build the RDataFrame
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df{*chain};

    // Gate on beam events
    auto gatedCFA = df.Filter([](ActRoot::ModularData &d){ return d.Get("GATCONF") == 64; }, {"ModularData"});
    auto gatedF2  = df.Filter([](ActRoot::ModularData &d){ return d.Get("GATCONF") == 32; }, {"ModularData"});

    // Encontrar rango de runs
    auto minRun = *df.Min<int>("fRun");
    auto maxRun = *df.Max<int>("fRun");
    int nbinsX = maxRun - minRun + 1;
    std::cout << "Rango de runs: " << minRun << " - " << maxRun 
              << " → nbinsX = " << nbinsX << std::endl;

    // Histogramas por run (un bin por run)
    auto cfa_counts = gatedCFA.Histo1D({"cfaCounts", "CFA counts per run;Run;Counts",
                                        nbinsX, minRun - 0.5, maxRun + 0.5}, "fRun");
    auto f2_counts  = gatedF2.Histo1D({"f2Counts", "F2 counts per run;Run;Counts",
                                        nbinsX, minRun - 0.5, maxRun + 0.5}, "fRun");

    // Clonar para calcular el cociente
    auto cfa_h = (TH1D*)cfa_counts->Clone("cfa_h");
    auto f2_h  = (TH1D*)f2_counts->Clone("f2_h");

    auto ratio_h = (TH1D*)cfa_h->Clone("ratio_h");
    ratio_h->SetTitle("CFA/F2 counts per run;Run;CFA/F2");

    // Calcular cociente bin a bin
    for (int i = 1; i <= ratio_h->GetNbinsX(); i++) {
        double f2_val = f2_h->GetBinContent(i);
        if (f2_val != 0) {
            ratio_h->SetBinContent(i, cfa_h->GetBinContent(i)/f2_val);
        } else {
            ratio_h->SetBinContent(i, 0); // o NaN si quieres marcarlo
        }
    }

    // Normalizar al primer run
    double firstBin = ratio_h->GetBinContent(1);
    if (firstBin != 0) {
        ratio_h->Scale(1.0 / firstBin);
    }
    // Crear TGraph a partir del histograma
    int nbins = ratio_h->GetNbinsX();
    TGraph *gRatio = new TGraph(nbins);
    for (int i = 1; i <= nbins; i++) {
        double x = ratio_h->GetBinCenter(i);
        double y = ratio_h->GetBinContent(i);
        gRatio->SetPoint(i-1, x, y);
    }

    

    // Canvas
    auto c1 = new TCanvas("c1", "Beam Counts per Run", 1200, 800);
    c1->Divide(1,2);

    // Pad 1: CFA y F2
    c1->cd(1);
    cfa_h->SetLineColor(kRed);
    cfa_h->SetMarkerColor(kRed);
    cfa_h->SetMarkerStyle(20);
    cfa_h->DrawClone(); // Draw con errores
    f2_h->SetLineColor(kBlue);
    f2_h->SetMarkerColor(kBlue);
    f2_h->SetMarkerStyle(21);
    f2_h->DrawClone("same"); // Misma pad
    auto leg = new TLegend(0.7,0.7,0.9,0.9);
    leg->AddEntry(cfa_h, "CFA", "lep");
    leg->AddEntry(f2_h, "F2", "lep");
    leg->Draw();

    // Pad 2: cociente CFA/F2
    c1->cd(2);
    gRatio->SetMarkerColor(kBlack);
    gRatio->SetMarkerStyle(22);
    gRatio->SetTitle("CFA/F2 counts per run (normalized to first run);Run;CFA/F2"); // título y ejes
    gRatio->Draw("AP"); // A = axis, P = points
}


