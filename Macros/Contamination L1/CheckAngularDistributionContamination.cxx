#include "ActMergerData.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"

#include <iostream>
#include <string>

// Compara la theta light (lab) de dos ficheros de entrada:
//   d: eventos con fLight.IsFilled() == false         (fichero tree_ex_<beam>_d_d_filtered.root)
//   p: eventos con Ex en [exMinP, exMaxP]             (fichero tree_ex_<beam>_d_p_filtered.root)
// normalize = true -> ambos histogramas normalizados a area 1 (para comparar formas)
// Canvas 1: overlay de theta + ratio d / p
// Canvas 2: cinematica theta Lab vs EVertex (d y p)
void CheckAngularDistributionContamination(bool normalize = true)
{
    std::string beam {"7Li"};

    // Seleccion en Ex para el p
    const double exMinP {5.};
    const double exMaxP {6.};

    // Binning de theta [deg]
    const int nTh {180};
    const double thMin {0.}, thMax {180.};

    // Binning de EVertex [MeV] para la cinematica (ajusta el rango a tu haz)
    const int nE {20};
    const double eMin {0.}, eMax {5.};

    // ---- Ficheros de entrada
    ROOT::RDataFrame dfD {"Final_Tree",
                          TString::Format("../../PostAnalysis/Outputs/tree_ex_%s_d_d_filtered.root", beam.c_str())};
    ROOT::RDataFrame dfP {"Final_Tree",
                          TString::Format("../../PostAnalysis/Outputs/tree_ex_%s_d_p_filtered.root", beam.c_str())};

    // ---- d: solo eventos con fLight no filled
    auto dfDsel = dfD.Filter([](const ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"})
                      .Define("thL", [](const ActRoot::MergerData& m) { return static_cast<double>(m.fThetaLight); },
                              {"MergerData"});

    // ---- p: solo Ex en [5, 6]
    auto dfPsel = dfP.Filter([](const ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"})
                      .Filter([=](double ex) { return ex >= exMinP && ex <= exMaxP; }, {"Ex"})
                      .Define("thL", [](const ActRoot::MergerData& m) { return static_cast<double>(m.fThetaLight); },
                              {"MergerData"});

    // ---- Histogramas (lazy)
    auto hDres = dfDsel.Histo1D({"hThetaD", ";#theta_{Light} [#circ];Counts", nTh, thMin, thMax}, "thL");
    auto hPres = dfPsel.Histo1D({"hThetaP", ";#theta_{Light} [#circ];Counts", nTh, thMin, thMax}, "thL");

    // ---- Cinematica: theta Lab vs EVertex (lazy)
    auto hD2res = dfDsel.Histo2D(
        {"hKinD", "d ;#theta_{Light} [#circ];E_{Vertex} [MeV]", nTh, thMin, thMax, nE, eMin, eMax}, "thL",
        "EVertex");
    auto hP2res = dfPsel.Histo2D(
        {"hKinP", "p;#theta_{Light} [#circ];E_{Vertex} [MeV]", nTh, thMin, thMax, nE, eMin, eMax}, "thL", "EVertex");

    // Un solo trigger para los dos event loops
    // ROOT::RDF::RunGraphs({hDres, hPres, hD2res, hP2res});

    auto* hD = static_cast<TH1D*>(hDres->Clone("hThetaD_clone"));
    auto* hP = static_cast<TH1D*>(hPres->Clone("hThetaP_clone"));
    hD->SetDirectory(nullptr);
    hP->SetDirectory(nullptr);

    const double nD {hD->Integral()};
    const double nP {hP->Integral()};
    std::cout << "d (IsFilled == false): " << nD << " eventos\n";
    std::cout << "p (" << exMinP << " < Ex < " << exMaxP << " MeV): " << nP << " eventos\n";

    if(normalize)
    {
        if(nD > 0)
            hD->Scale(1. / nD);
        if(nP > 0)
            hP->Scale(1. / nP);
        hD->GetYaxis()->SetTitle("Normalized counts");
    }

    hD->SetLineColor(kAzure + 1);
    hD->SetLineWidth(2);
    hP->SetLineColor(kRed + 1);
    hP->SetLineWidth(2);

    // ---- Canvas 1: overlay + ratio d / p
    gStyle->SetOptStat(0);
    auto* c = new TCanvas("cThetaLight", "Theta light: d vs p", 900, 800);

    auto* pTop = new TPad("pTop", "pTop", 0, 0.3, 1, 1);
    pTop->SetBottomMargin(0.02);
    pTop->Draw();
    auto* pBot = new TPad("pBot", "pBot", 0, 0, 1, 0.3);
    pBot->SetTopMargin(0.03);
    pBot->SetBottomMargin(0.3);
    pBot->Draw();

    pTop->cd();
    hD->GetXaxis()->SetLabelSize(0);
    hD->GetXaxis()->SetTitleSize(0);
    hD->SetMaximum(1.15 * TMath::Max(hD->GetMaximum(), hP->GetMaximum()));
    hD->SetMinimum(0);
    hD->DrawClone("hist");
    hP->DrawClone("hist same");

    auto* leg = new TLegend(0.60, 0.72, 0.89, 0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(hD, TString::Format("d (not filled), N = %.0f", nD), "l");
    leg->AddEntry(hP, TString::Format("p (%.0f < E_{x} < %.0f MeV), N = %.0f", exMinP, exMaxP, nP), "l");
    leg->Draw();

    pBot->cd();
    auto* hR = static_cast<TH1D*>(hD->Clone("hRatio"));
    hR->SetDirectory(nullptr);
    hR->Divide(hP);
    hR->SetTitle("");
    hR->GetYaxis()->SetTitle("d / p");
    hR->GetYaxis()->SetNdivisions(505);
    hR->GetYaxis()->SetTitleSize(0.10);
    hR->GetYaxis()->SetTitleOffset(0.45);
    hR->GetYaxis()->SetLabelSize(0.09);
    hR->GetXaxis()->SetTitle("#theta_{Light} [#circ]");
    hR->GetXaxis()->SetLabelSize(0.10);
    hR->GetXaxis()->SetTitleSize(0.12);
    hR->GetXaxis()->SetTitleOffset(1.0);
    hR->SetLineColor(kBlack);
    hR->SetMarkerStyle(20);
    hR->SetMarkerSize(0.6);
    hR->SetMinimum(0);
    hR->SetMaximum(3);
    hR->DrawClone("pe");

    // ---- Canvas 2: cinematica theta Lab vs EVertex
    auto* hKinD = static_cast<TH2D*>(hD2res->Clone("hKinD_clone"));
    auto* hKinP = static_cast<TH2D*>(hP2res->Clone("hKinP_clone"));
    hKinD->SetDirectory(nullptr);
    hKinP->SetDirectory(nullptr);

    auto* cK = new TCanvas("cKinematics", "Kinematics: theta Light vs EVertex", 1400, 600);
    cK->Divide(2, 1);

    cK->cd(1);
    gPad->SetRightMargin(0.14);
    gPad->SetLogz();
    hKinD->SetTitle(TString::Format("d (not filled), N = %.0f", nD));
    hKinD->Draw("colz");

    cK->cd(2);
    gPad->SetRightMargin(0.14);
    gPad->SetLogz();
    hKinP->SetTitle(TString::Format("p (%.0f < E_{x} < %.0f MeV), N = %.0f", exMinP, exMaxP, nP));
    hKinP->Draw("colz");
}