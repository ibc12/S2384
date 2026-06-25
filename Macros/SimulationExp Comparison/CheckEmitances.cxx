#include "ActInputParser.h"
#include "ActMergerData.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TString.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>


// Macro to check all the 11Li emmitances, for all the beam positions
void CheckEmitances()
{
    // Get all the emmitances from the simulation Input folder
    std::string emittanceFirst = "../Emittance/Outputs/emittance11Li_pre.root";
    std::string emittanceSecond = "../Emittance/Outputs/emittance11Li_post_preMeasure.root";
    std::string emittanceThird = "../Emittance/Outputs/emittance11Li_post_postMeasure.root";

    // Get the .root files and check the histograms
    TFile* f1 = TFile::Open(emittanceFirst.c_str());
    TFile* f2 = TFile::Open(emittanceSecond.c_str());
    TFile* f3 = TFile::Open(emittanceThird.c_str());

    // Get dataframe Emittance_Tree
    ROOT::RDataFrame df1("Emittance_Tree", f1);
    ROOT::RDataFrame df2("Emittance_Tree", f2);
    ROOT::RDataFrame df3("Emittance_Tree", f3);

    // Get histos of AtBegin.fCoordinates.fY with bining from 110 to 150 with 40 bins
    auto hY1begin = df1.Histo1D({"hY1begin", "Y at Begin;Y [mm];Counts", 160, 110, 150}, "AtBegin.fCoordinates.fY");
    auto hY2begin = df2.Histo1D({"hY2begin", "Y at Begin;Y [mm];Counts", 160, 110, 150}, "AtBegin.fCoordinates.fY");
    auto hY3begin = df3.Histo1D({"hY3begin", "Y at Begin;Y [mm];Counts", 160, 110, 150}, "AtBegin.fCoordinates.fY");
    auto hY1end = df1.Histo1D({"hY1end", "Y at End;Y [mm];Counts", 160, 110, 150}, "AtEnd.fCoordinates.fY");
    auto hY2end = df2.Histo1D({"hY2end", "Y at End;Y [mm];Counts", 160, 110, 150}, "AtEnd.fCoordinates.fY");
    auto hY3end = df3.Histo1D({"hY3end", "Y at End;Y [mm];Counts", 160, 110, 150}, "AtEnd.fCoordinates.fY");

    // Fit histos to gaussian
    hY1begin->Fit("gaus");
    hY2begin->Fit("gaus");
    hY3begin->Fit("gaus");
    hY1end->Fit("gaus");
    hY2end->Fit("gaus");
    hY3end->Fit("gaus");

    // Draw the histograms
    TCanvas* c1 = new TCanvas("c1", "Emittance Comparison", 800, 600);
    hY1begin->SetLineColor(kRed);
    hY2begin->SetLineColor(kBlue);
    hY3begin->SetLineColor(kGreen);
    hY1begin->DrawClone();
    hY2begin->DrawClone("SAME");
    hY3begin->DrawClone("SAME");
    // Draw result of fit
    hY1begin->GetFunction("gaus")->SetLineColor(kRed);
    hY2begin->GetFunction("gaus")->SetLineColor(kBlue);
    hY3begin->GetFunction("gaus")->SetLineColor(kGreen);
    // Write result of fit in canvas
    TLatex* latex1 = new TLatex(0.6, 0.8,
                                Form("Mean: %.2f, Sigma: %.2f", hY1begin->GetFunction("gaus")->GetParameter(1),
                                     hY1begin->GetFunction("gaus")->GetParameter(2)));
    latex1->SetNDC();
    latex1->SetTextColor(kRed);
    latex1->Draw();
    TLatex* latex2 = new TLatex(0.6, 0.75,
                                Form("Mean: %.2f, Sigma: %.2f", hY2begin->GetFunction("gaus")->GetParameter(1),
                                     hY2begin->GetFunction("gaus")->GetParameter(2)));
    latex2->SetNDC();
    latex2->SetTextColor(kBlue);
    latex2->Draw();
    TLatex* latex3 = new TLatex(0.6, 0.7,
                                Form("Mean: %.2f, Sigma: %.2f", hY3begin->GetFunction("gaus")->GetParameter(1),
                                     hY3begin->GetFunction("gaus")->GetParameter(2)));
    latex3->SetNDC();
    latex3->SetTextColor(kGreen);
    latex3->Draw();

    // Draw for the end histograms
    TCanvas* c2 = new TCanvas("c2", "Emittance Comparison End", 800, 600);
    hY1end->SetLineColor(kRed);
    hY2end->SetLineColor(kBlue);
    hY3end->SetLineColor(kGreen);
    hY1end->DrawClone();
    hY2end->DrawClone("SAME");
    hY3end->DrawClone("SAME");
    // Draw result of fit
    hY1end->GetFunction("gaus")->SetLineColor(kRed);
    hY2end->GetFunction("gaus")->SetLineColor(kBlue);
    hY3end->GetFunction("gaus")->SetLineColor(kGreen);
    // Write result of fit in canvas
    TLatex* latex4 = new TLatex(0.6, 0.8,
                                Form("Mean: %.2f, Sigma: %.2f", hY1end->GetFunction("gaus")->GetParameter(1),
                                     hY1end->GetFunction("gaus")->GetParameter(2)));
    latex4->SetNDC();
    latex4->SetTextColor(kRed);
    latex4->Draw();
    TLatex* latex5 = new TLatex(0.6, 0.75,
                                Form("Mean: %.2f, Sigma: %.2f", hY2end->GetFunction("gaus")->GetParameter(1),
                                     hY2end->GetFunction("gaus")->GetParameter(2)));
    latex5->SetNDC();
    latex5->SetTextColor(kBlue);
    latex5->Draw();
    TLatex* latex6 = new TLatex(0.6, 0.7,
                                Form("Mean: %.2f, Sigma: %.2f", hY3end->GetFunction("gaus")->GetParameter(1),
                                     hY3end->GetFunction("gaus")->GetParameter(2)));
    latex6->SetNDC();
    latex6->SetTextColor(kGreen);
    latex6->Draw();
}