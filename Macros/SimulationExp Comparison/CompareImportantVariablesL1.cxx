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
#include <cmath>
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

// Reconstruct RP.X from the end of the light track (beam along X):
// RP.X = LastPosX - TL * cos(theta_lab). LastPosX and TL in mm, theta in deg.
inline double ReconstructRPx(double lastX, double tl, double thetaDeg)
{
    return lastX - tl * std::cos(thetaDeg * M_PI / 180.);
}

// Insert pads (unique X, Y) outside the exclusion zone
template <typename Container>
void FillPadsOutZone(const Container& voxels, std::set<std::pair<int, int>>& pads)
{
    for(const auto& voxel : voxels)
    {
        auto pos = voxel.GetPosition();
        if(pos.Y() < yMinExclusionZone || pos.Y() > yMaxExclusionZone)
            pads.insert({static_cast<int>(pos.X()), static_cast<int>(pos.Y())});
    }
}

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

// With this macro we saw that the experiment rp.x ends about 240 mm while the simulated ends in 256. We have to adjust
// the simulation to match the experiment. But we also should

void CompareImportantVariablesL1()
{
    // Reaction: 7Li(d,p)8Li, only L1 events, elastic-like Ex = 0
    std::string beam {"7Li"};
    std::string target {"2H"};
    std::string targetExp {"d"};
    std::string light {"1H"};
    std::string lightExp {"p"};
    double Ex {0.}; // MeV

    // For the 2D RP.X vs theta plot of the experiment:
    //  false -> use the real RP.X (MergerData::fRP.X())
    //  true  -> use the same reconstruction as in the simulation (most consistent comparison)
    const bool useRecoRPExp {false};

    // TString fileNameSimu {TString::Format("../../Simulation/Outputs/%s/test_rpCut/"
    //                                       "%s_%s_TRIUMF_Eex_%.3f_nPS_%d_pPS_%d_L1_232rpCut.root",
    //                                       beam.c_str(), target.c_str(), light.c_str(), Ex, 0, 0)};
    TString fileNameSimu {TString::Format("../../Simulation/Outputs/%s/test_rpCut/"
                                          "%s_%s_TRIUMF_Eex_%.3f_nPS_%d_pPS_%d_L1_240rpCut_paralelization.root",
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

    // Simu columns used: theta3Lab, phi3CM, TL, RawTL, EVertex, nPads, LastPosX, LastPosY
    // RPx is reconstructed from the end of the light track (it is not a branch of the simulation tree)
    auto dfSimu = dfSimuRaw.Define("phiLabDeg", [](double phi) { return phi * 180.0 / M_PI; }, {"phi3CM"})
                      .Define("RPx", ReconstructRPx, {"LastPosX", "TL", "theta3Lab"});

    // ---------------- Experiment: only L1 events and unify column names ----------------
    auto dfExp =
        dfExpRaw.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"})
            .Filter([](double ex) { return ex < 0.5; }, {"Ex"})
            .Define("thetaLabDeg", [](ActRoot::MergerData& m) { return static_cast<double>(m.fThetaLight); },
                    {"MergerData"})
            .Define("phiLabDeg", [](ActRoot::MergerData& m) { return static_cast<double>(m.fPhiLight); },
                    {"MergerData"})
            .Define("RPx", [](ActRoot::MergerData& m) { return static_cast<double>(m.fRP.X()); }, {"MergerData"})
            .Define("TLExp", [](ActRoot::MergerData& m) { return static_cast<double>(m.fLight.fTL); }, {"MergerData"})
            .Define("RawTLExp", [](ActRoot::MergerData& m) { return static_cast<double>(m.fLight.fRawTL); },
                    {"MergerData"})
            // Last voxel: pads -> mm, to match the simulation
            .Define("LastPosX", [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                    { return padSize * GetLastVoxelXY(m, tpc).first; }, {"MergerData", "TPCData"})
            .Define("LastPosY", [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                    { return padSize * GetLastVoxelXY(m, tpc).second; }, {"MergerData", "TPCData"})
            // RP.X reconstructed with exactly the same formula as in the simulation
            .Define("RPxRec", ReconstructRPx, {"LastPosX", "TLExp", "thetaLabDeg"})
            .Define("nPadsTotal",
                    [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                    {
                        std::set<std::pair<int, int>> pads;
                        FillPadsOutZone(tpc.fClusters[m.fLightIdx].GetVoxels(), pads);
                        FillPadsOutZone(tpc.fClusters[m.fHeavyIdx].GetVoxels(), pads);
                        FillPadsOutZone(tpc.fClusters[m.fBeamIdx].GetVoxels(), pads);
                        FillPadsOutZone(tpc.fRaw, pads);
                        return static_cast<int>(pads.size());
                    },
                    {"MergerData", "TPCData"});

    // ---------------- Histograms ----------------
    // Total number of pads out of the exclusion zone
    auto hPadsSimu =
        dfSimu.Histo1D({"hPadsSimu", "Total pads out of exclusion zone;N pads;Counts", 150, 0, 150}, "nPads");
    auto hPadsExp = dfExp.Histo1D({"hPadsExp", "", 150, 0, 150}, "nPadsTotal");

    // Angular distributions
    auto hThetaSimu = dfSimu.Histo1D(
        {"hThetaSimu", "Theta lab;#theta_{Lab} (deg);Counts", thetaInterval, thetaMin, thetaMax}, "theta3Lab");
    auto hThetaExp = dfExp.Histo1D({"hThetaExp", "", thetaInterval, thetaMin, thetaMax}, "thetaLabDeg");
    auto hPhiSimu = dfSimu.Histo1D({"hPhiSimu", "Phi;#phi (deg);Counts", 72, -180, 180}, "phiLabDeg");
    auto hPhiExp = dfExp.Histo1D({"hPhiExp", "", 72, -180, 180}, "phiLabDeg");

    // Reaction point X: exp histogram uses the same binning as the simu hRP
    const int nBinsRP = hRPSimu->GetNbinsX();
    const double rpMin = hRPSimu->GetXaxis()->GetXmin();
    const double rpMax = hRPSimu->GetXaxis()->GetXmax();
    ROOT::RDF::TH1DModel rpModel {"hRPxExp", "", nBinsRP, rpMin, rpMax};
    auto hRPxExp = dfExp.Histo1D(rpModel, "RPx");

    // RP.X vs theta lab (simu: reconstructed RP.X; exp: real or reconstructed depending on useRecoRPExp)
    auto hRPThetaSimu = dfSimu.Histo2D({"hRPThetaSimu", "RP.X vs #theta_{Lab} (simu);#theta_{Lab} (deg);RP.X (mm)",
                                        thetaInterval, 0, thetaMax, nBinsRP, rpMin, rpMax},
                                       "theta3Lab", "RPx");
    auto hRPThetaExp = dfExp.Histo2D({"hRPThetaExp", "RP.X vs #theta_{Lab} (exp);#theta_{Lab} (deg);RP.X (mm)",
                                      thetaInterval, 0, thetaMax, nBinsRP, rpMin, rpMax},
                                     "thetaLabDeg", useRecoRPExp ? "RPxRec" : "RPx");

    auto hLastXThetaSimu = dfSimu.Histo2D(
        {"hLastXThetaSimu", ";#theta_{Lab} (deg);LastPosX (mm)", thetaInterval, 0, thetaMax, 128, 0, detMax},
        "theta3Lab", "LastPosX");
    auto hLastXThetaExp = dfExp.Histo2D(
        {"hLastXThetaExp", ";#theta_{Lab} (deg);LastPosX (mm)", thetaInterval, 0, thetaMax, 128, 0, detMax},
        "thetaLabDeg", "LastPosX");

    // TL and RawTL
    auto hTLSimu = dfSimu.Histo1D({"hTLSimu", "TL;TL (mm);Counts", 100, 0, 200}, "TL");
    auto hTLExp = dfExp.Histo1D({"hTLExp", "", 100, 0, 200}, "TLExp");
    auto hRawTLSimu = dfSimu.Histo1D({"hRawTLSimu", "RawTL;RawTL (mm);Counts", 100, 0, 200}, "RawTL");
    auto hRawTLExp = dfExp.Histo1D({"hRawTLExp", "", 100, 0, 200}, "RawTLExp");

    // Last voxel position
    auto hLastSimu =
        dfSimu.Histo2D({"hLastSimu", "Last voxel (simu);LastPosX (mm);LastPosY (mm)", 128, 0, detMax, 128, 0, detMax},
                       "LastPosX", "LastPosY");
    auto hLastExp =
        dfExp.Histo2D({"hLastExp", "Last voxel (exp);LastPosX (mm);LastPosY (mm)", 128, 0, detMax, 128, 0, detMax},
                      "LastPosX", "LastPosY");

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
    DrawOverlay(hPadsSimu.GetPtr(), hPadsExp.GetPtr());
    c1->cd(2);
    DrawOverlay(hRPSimu, hRPxExp.GetPtr());
    c1->cd(3);
    DrawOverlay(hThetaSimu.GetPtr(), hThetaExp.GetPtr());
    c1->cd(4);
    DrawOverlay(hPhiSimu.GetPtr(), hPhiExp.GetPtr());

    auto c2 = new TCanvas("c2", "TL, RawTL and last voxel", 1200, 800);
    c2->Divide(2, 2);
    c2->cd(1);
    DrawOverlay(hTLSimu.GetPtr(), hTLExp.GetPtr());
    c2->cd(2);
    DrawOverlay(hRawTLSimu.GetPtr(), hRawTLExp.GetPtr());
    c2->cd(3);
    hLastSimu->DrawClone("colz");
    c2->cd(4);
    hLastExp->DrawClone("colz");

    auto c3 = new TCanvas("c3", "Kinematics", 1200, 500);
    c3->Divide(2, 1);
    c3->cd(1);
    hKinSimu->DrawClone("colz");
    c3->cd(2);
    hKinExp->DrawClone("colz");

    // RP.X vs theta lab: both normalised to unit area and sharing the same colour scale
    auto c4 = new TCanvas("c4", "RP.X and LastPoint vs Theta", 1200, 500);
    c4->Divide(2, 2);
    auto* hSimu2D = hRPThetaSimu.GetPtr();
    auto* hExp2D = hRPThetaExp.GetPtr();
    hSimu2D->SetStats(false);
    hExp2D->SetStats(false);
    if(hSimu2D->Integral() > 0)
        hSimu2D->Scale(1. / hSimu2D->Integral());
    if(hExp2D->Integral() > 0)
        hExp2D->Scale(1. / hExp2D->Integral());
    const double zMax = std::max(hSimu2D->GetMaximum(), hExp2D->GetMaximum());
    hSimu2D->SetMaximum(zMax);
    hExp2D->SetMaximum(zMax);
    c4->cd(1);
    hSimu2D->DrawClone("colz");
    c4->cd(2);
    hExp2D->DrawClone("colz");
    auto* hSimu2DLast = hLastXThetaSimu.GetPtr();
    auto* hExp2DLast = hLastXThetaExp.GetPtr();
    hSimu2DLast->SetStats(false);
    hExp2DLast->SetStats(false);
    if(hSimu2DLast->Integral() > 0)
        hSimu2DLast->Scale(1. / hSimu2DLast->Integral());
    if(hExp2DLast->Integral() > 0)
        hExp2DLast->Scale(1. / hExp2DLast->Integral());
    const double zMaxLast = std::max(hSimu2DLast->GetMaximum(), hExp2DLast->GetMaximum());
    hSimu2DLast->SetMaximum(zMaxLast);
    hExp2DLast->SetMaximum(zMaxLast);
    c4->cd(3);
    hSimu2DLast->DrawClone("colz");
    c4->cd(4);
    hExp2DLast->DrawClone("colz");
}