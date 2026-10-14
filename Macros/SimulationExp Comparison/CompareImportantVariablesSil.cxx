#include "ActMergerData.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"

#include <algorithm>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <utility>

// Exclusion zone in Y (pads)
constexpr int yMinExclusionZone = 56;
constexpr int yMaxExclusionZone = 71;

// Histogram ranges (adjust if needed)
constexpr double padSize = 2.;                     // mm per pad
constexpr double detMax = 256.;                    // mm (last voxel histograms)
constexpr double thetaMin = 130.;                  // deg
constexpr double thetaMax = 170.;                  // deg
constexpr int thetaInterval = thetaMax - thetaMin; // deg
constexpr double eVertexMax = 3.;                  // MeV

// Position (X, Y) in pads of the last voxel of the light cluster
std::pair<double, double> GetLastVoxelXY(ActRoot::MergerData& m, ActRoot::TPCData& tpc)
{
    auto lightCl = tpc.fClusters[m.fLightIdx]; // copy: we modify it
    auto line = lightCl.GetLine();
    auto vertex = tpc.fRPs.front();
    auto directionToSort = line.GetPoint() - vertex;
    // Sort voxels according to their projection on the line
    lightCl.SortAlongDir(directionToSort);
    lightCl.ReFit();
    auto pos = lightCl.GetVoxels().back().GetPosition();
    return {pos.X(), pos.Y()};
}

// Overlay simu (red, line) and exp (black, points), both normalised to unit area
void DrawOverlay(TH1* hSimu, TH1* hExp)
{
    hSimu->SetStats(false);
    hExp->SetStats(false);
    hSimu->SetLineColor(kRed);
    hSimu->SetLineWidth(2);
    hExp->SetLineColor(kBlack);
    hExp->SetMarkerColor(kBlack);
    hExp->SetMarkerStyle(20);
    hExp->SetMarkerSize(0.6);

    if(hSimu->Integral() > 0)
        hSimu->Scale(1. / hSimu->Integral());
    if(hExp->Integral() > 0)
        hExp->Scale(1. / hExp->Integral());

    hSimu->SetMaximum(1.2 * std::max(hSimu->GetMaximum(), hExp->GetMaximum()));
    hSimu->GetYaxis()->SetTitle("Normalised counts");

    hSimu->DrawClone("HIST");
    hExp->DrawClone("E SAME");

    auto leg = new TLegend(0.65, 0.75, 0.88, 0.88);
    leg->AddEntry(hSimu, "Simu", "l");
    leg->AddEntry(hExp, "Exp", "lep");
    leg->Draw();
}

void CompareImportantVariablesSil()
{
    // Reaction: 7Li(d,p)8Li, only L1 events, elastic-like Ex = 0
    std::string beam {"7Li"};
    std::string target {"2H"};
    std::string targetExp {"d"};
    std::string light {"1H"};
    std::string lightExp {"p"};
    double Ex {0.}; // MeV

    TString fileNameSimu {TString::Format("../../Simulation/Outputs/%s/"
                                          "%s_%s_TRIUMF_Eex_%.3f_nPS_%d_pPS_%d.root",
                                          beam.c_str(), target.c_str(), light.c_str(), Ex, 0, 0)};
    TString fileNameExp {TString::Format("../../PostAnalysis/Outputs/tree_ex_F_%s_%s_%s.root", beam.c_str(),
                                         targetExp.c_str(), lightExp.c_str())};

    ROOT::EnableImplicitMT();
    gStyle->SetOptStat(0);

    // ---------------- Simulation: RP is not a branch, it is saved as histogram "hRP" in the TFile ----------------
    std::unique_ptr<TFile> fSimu {TFile::Open(fileNameSimu)};
    if(!fSimu || fSimu->IsZombie())
    {
        std::cout << "Could not open simulation file " << fileNameSimu << '\n';
        return;
    }
    // hRP is a TH2: RP.X() is on the X axis and RP.Y() on the Y axis (whatever the axis titles say)
    auto* hRP2D = fSimu->Get<TH2>("hRP");
    if(!hRP2D)
    {
        std::cout << "Histogram hRP (TH2) not found in " << fileNameSimu << '\n';
        return;
    }
    // Project onto X (RP.X) and detach from the file so it survives after closing it
    auto* hRPSimu = hRP2D->ProjectionX("hRPSimu");
    hRPSimu->SetDirectory(nullptr);
    hRPSimu->SetTitle("RP.X;RP.X;Counts");
    fSimu->Close();

    auto dfSimuRaw = ROOT::RDataFrame("SimulationTTree", fileNameSimu);
    auto dfExpRaw = ROOT::RDataFrame("Final_Tree", fileNameExp);

    // Simu columns used: theta3Lab, phi3CM, TL, RawTL, EVertex, LastPosX, LastPosY
    auto dfSimu = dfSimuRaw.Define("phiLabDeg", [](double phi) { return phi * 180.0 / M_PI; }, {"phi3CM"});

    // ---------------- Experiment: only L1 events and unify column names ----------------
    auto dfExp =
        dfExpRaw.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == true; }, {"MergerData"})
            .Filter([](double ex) { return ex < 0.5; }, {"Ex"})
            .Define("thetaLabDeg", [](ActRoot::MergerData& m) { return static_cast<double>(m.fThetaLight); },
                    {"MergerData"})
            .Define("phiLabDeg", [](ActRoot::MergerData& m) { return static_cast<double>(m.fPhiLight); },
                    {"MergerData"})
            .Define("RPx", [](ActRoot::MergerData& m) { return static_cast<double>(m.fRP.X()); }, {"MergerData"})
            .Define("TLExp", [](ActRoot::MergerData& m) { return static_cast<double>(m.fLight.fTL); }, {"MergerData"})
            .Define("RawTLExp", [](ActRoot::MergerData& m) { return static_cast<double>(m.fLight.fRawTL); },
                    {"MergerData"});

    // ---------------- Histograms ----------------

    // Angular distributions
    auto hThetaSimu = dfSimu.Histo1D(
        {"hThetaSimu", "Theta lab;#theta_{Lab} (deg);Counts", thetaInterval, thetaMin, thetaMax}, "theta3Lab");
    auto hThetaExp = dfExp.Histo1D({"hThetaExp", "", thetaInterval, thetaMin, thetaMax}, "thetaLabDeg");
    auto hPhiSimu = dfSimu.Histo1D({"hPhiSimu", "Phi;#phi (deg);Counts", 72, -180, 180}, "phiLabDeg");
    auto hPhiExp = dfExp.Histo1D({"hPhiExp", "", 72, -180, 180}, "phiLabDeg");

    // Reaction point X: exp histogram uses the same binning as the simu hRP
    ROOT::RDF::TH1DModel rpModel {"hRPxExp", "", hRPSimu->GetNbinsX(), hRPSimu->GetXaxis()->GetXmin(),
                                  hRPSimu->GetXaxis()->GetXmax()};
    auto hRPxExp = dfExp.Histo1D(rpModel, "RPx");

    // Kinematics: thetaLab vs EVertex
    auto hKinSimu = dfSimu.Histo2D({"hKinSimu", "Kinematics (simu);#theta_{Lab} (deg);E_{Vertex} (MeV)", thetaInterval,
                                    thetaMin, thetaMax, 120, 0, eVertexMax},
                                   "theta3Lab", "EVertex");
    auto hKinExp = dfExp.Histo2D({"hKinExp", "Kinematics (exp);#theta_{Lab} (deg);E_{Vertex} (MeV)", thetaInterval,
                                  thetaMin, thetaMax, 120, 0, eVertexMax},
                                 "thetaLabDeg", "EVertex");

    // ---------------- Plots ----------------
    auto c1 = new TCanvas("c1", "Pads, angles and RP.X", 1200, 800);
    c1->Divide(2, 2);
    c1->cd(1);
    c1->cd(2);
    DrawOverlay(hRPSimu, hRPxExp.GetPtr());
    c1->cd(3);
    DrawOverlay(hThetaSimu.GetPtr(), hThetaExp.GetPtr());
    c1->cd(4);
    DrawOverlay(hPhiSimu.GetPtr(), hPhiExp.GetPtr());

    auto c2 = new TCanvas("c2", "TL, RawTL and last voxel", 1200, 800);
    c2->Divide(2, 2);
    c2->cd(1);
    c2->cd(3);

    auto c3 = new TCanvas("c3", "Kinematics", 1200, 500);
    c3->Divide(2, 1);
    c3->cd(1);
    hKinSimu->DrawClone("colz");
    c3->cd(2);
    hKinExp->DrawClone("colz");
}