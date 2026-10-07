#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActParticle.h"
#include "ActRunner.h"
#include "ActSRIM.h"
#include "ActSilData.h"
#include "ActSilMatrix.h"
#include "ActSilSpecs.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TH2.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TMath.h"
#include "TObject.h"
#include "TRandom.h"
#include "TString.h"

#include "Math/AxisAngle.h"
#include "Math/DisplacementVector3D.h"
#include "Math/Point3Dfwd.h"
#include "Math/Rotation3D.h"
#include "Math/RotationZYX.h"
#include "Math/Vector3D.h"

#include <fstream>
#include <map>
#include <string>

#include "../../PostAnalysis/HistConfig.h"

using XYZPoint = ROOT::Math::XYZPointF;
using XYZVector = ROOT::Math::XYZVectorF;

// double GetPhi3D(const XYZVector& beam, const XYZVector& other)
// {
//     auto beamFrame {beam.Unit()};
//     XYZVector worldFrame {1, 0, 0};
//
//     // Rotación inversa a la de la simulación: lleva el haz -> X (lab -> frame del haz)
//     // Eje invertido (u x X en vez de X x u), mismo ángulo
//     auto cross {beamFrame.Cross(worldFrame)};
//     auto angle {TMath::ACos(beamFrame.Dot(worldFrame))};
//     ROOT::Math::AxisAngle axis {cross, angle};
//     ROOT::Math::Rotation3D rotation {axis};
//
//     // Traza en el frame del haz
//     auto t {rotation(other.Unit())};
//     return TMath::ATan2(t.Y(), t.Z()) * TMath::RadToDeg();
// }

double GetPhi3D(const XYZVector& beam, const XYZVector& other)
{
    const auto b {beam.Unit()};
    const XYZVector ex {1, 0, 0};

    const auto cross {b.Cross(ex)};
    const double s {cross.R()}; // sin(angulo)
    const double c {b.Dot(ex)}; // cos(angulo)

    ROOT::Math::Rotation3D rot; // identidad por defecto
    if(s > 1e-12)
        rot = ROOT::Math::Rotation3D {ROOT::Math::AxisAngle {cross / s, std::atan2(s, c)}};
    else if(c < 0) // haz antiparalelo a X: giro de pi sobre un eje perpendicular
        rot = ROOT::Math::Rotation3D {ROOT::Math::AxisAngle {XYZVector {0, 0, 1}, TMath::Pi()}};

    const auto t {rot(other.Unit())};
    return TMath::ATan2(t.Y(), t.Z()) * TMath::RadToDeg();
}


double GetPhi3DLegacy(const XYZVector& beam, const XYZVector& other)
{
    // TODO: Check validity of phi calculation

    // auto ub {beam.Unit()};            // unitary beam
    auto trackUnitary {other.Unit()};
    // XYZVector yz {0, ub.Y(), ub.Z()}; // beam dir in YZ plane
    // auto dot {other.Unit().Dot(yz) / yz.R()};
    // return TMath::ACos(dot) * TMath::RadToDeg();
    return TMath::ATan2(trackUnitary.Y(), trackUnitary.Z()) * TMath::RadToDeg();
}

void CheckPhiSumDistribution()
{
    std::string beam {"7Li"};
    std::string light {"p"};
    // Get df from root file
    ROOT::RDataFrame dfini {"Final_Tree", TString::Format("../../PostAnalysis/Outputs/tree_ex_%s_d_%s_filtered.root",
                                                          beam.c_str(), light.c_str())};

    double Sn_7Li {7.};   // MeV
    double Sn_8Li {2.03}; // MeV

    // Get drift parameter
    ActRoot::InputParser parser {};
    parser.ReadFile("../../configs/detector.conf");
    auto driftBlock = parser.GetBlock("Merger");
    auto driftFactor = driftBlock->GetDouble("DriftFactor");

    auto df = dfini
                  .Define("phiBeam",
                          [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                          {
                              // Get phi of beam for reference
                              int beamIdx = m.fBeamIdx;
                              auto dir = tpc.fClusters[beamIdx].GetLine().GetDirection().Unit();
                              return TMath::ATan2(dir.Y(), dir.Z()) *
                                     TMath::RadToDeg(); // Same definition as in ActMergerDetector
                          },
                          {"MergerData", "TPCData"})
                  .Define("phiSum",
                          [](ActRoot::MergerData& m, double phiBeam)
                          {
                              double substraction {0};
                              if(m.fPhiLight > m.fPhiHeavy)
                                  substraction += (m.fPhiLight - m.fPhiHeavy);
                              else
                                  substraction += (-m.fPhiHeavy + m.fPhiLight);
                              return substraction - 180;
                          },
                          {"MergerData", "phiBeam"})
                  .Define("phiLightLegacy",
                          [&driftFactor](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                          {
                              int lightIdx = m.fLightIdx;
                              int beamIdx = m.fBeamIdx;
                              auto line = tpc.fClusters[lightIdx].GetLine();
                              auto beamLine = tpc.fClusters[beamIdx].GetLine();
                              line.Scale(2, driftFactor);
                              beamLine.Scale(2, driftFactor);
                              auto dir = line.GetDirection();
                              auto beamDir = beamLine.GetDirection();
                              return GetPhi3DLegacy(beamDir, dir);
                          },
                          {"MergerData", "TPCData"})
                  .Define("phiHeavyLegacy",
                          [&driftFactor](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                          {
                              int heavyIdx = m.fHeavyIdx;
                              int beamIdx = m.fBeamIdx;
                              auto line = tpc.fClusters[heavyIdx].GetLine();
                              auto beamLine = tpc.fClusters[beamIdx].GetLine();
                              line.Scale(2, driftFactor);
                              beamLine.Scale(2, driftFactor);
                              auto dir = line.GetDirection();
                              auto beamDir = beamLine.GetDirection();
                              return GetPhi3DLegacy(beamDir, dir);
                          },
                          {"MergerData", "TPCData"})
                  .Define("phiLight",
                          [&driftFactor](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                          {
                              int lightIdx = m.fLightIdx;
                              int beamIdx = m.fBeamIdx;
                              auto line = tpc.fClusters[lightIdx].GetLine();
                              auto beamLine = tpc.fClusters[beamIdx].GetLine();
                              line.Scale(2, driftFactor);
                              beamLine.Scale(2, driftFactor);
                              auto dir = line.GetDirection();
                              auto beamDir = beamLine.GetDirection();
                              return GetPhi3D(beamDir, dir);
                          },
                          {"MergerData", "TPCData"})
                  .Define("phiHeavy",
                          [&driftFactor](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                          {
                              int heavyIdx = m.fHeavyIdx;
                              int beamIdx = m.fBeamIdx;
                              auto line = tpc.fClusters[heavyIdx].GetLine();
                              auto beamLine = tpc.fClusters[beamIdx].GetLine();
                              line.Scale(2, driftFactor);
                              beamLine.Scale(2, driftFactor);
                              auto dir = line.GetDirection();
                              auto beamDir = beamLine.GetDirection();
                              return GetPhi3D(beamDir, dir);
                          },
                          {"MergerData", "TPCData"})
                  .Define("phiSumLegacy",
                          [](double phiLightLegacy, double phiHeavyLegacy)
                          {
                              double substraction {0};
                              if(phiLightLegacy > phiHeavyLegacy)
                                  substraction += (phiLightLegacy - phiHeavyLegacy);
                              else
                                  substraction += (-phiHeavyLegacy + phiLightLegacy);
                              return substraction - 180;
                          },
                          {"phiLightLegacy", "phiHeavyLegacy"})
                  .Define("phiSumNew",
                          [](double phiLight, double phiHeavy)
                          {
                              double substraction {0};
                              if(phiLight > phiHeavy)
                                  substraction += (phiLight - phiHeavy);
                              else
                                  substraction += (-phiHeavy + phiLight);
                              return substraction - 180;
                          },
                          {"phiLight", "phiHeavy"});

    auto df_sil = df.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == true; }, {"MergerData"});
    auto df_sil_bound = df_sil.Filter(
        [&light, &Sn_7Li, &Sn_8Li](double ex)
        {
            if(light == "d" && ex < Sn_7Li)
                return true;
            else if(light == "p" && ex < Sn_8Li)
                return true;
            else
                return false;
        },
        {"Ex"});
    auto df_sil_unbound = df_sil.Filter(
        [&light, &Sn_7Li, &Sn_8Li](double ex)
        {
            if(light == "d" && ex > Sn_7Li)
                return true;
            else if(light == "p" && ex > Sn_8Li)
                return true;
            else
                return false;
        },
        {"Ex"});
    auto df_L1 = df.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"});
    auto df_L1_bound = df_L1.Filter(
        [&light, &Sn_7Li, &Sn_8Li](double ex)
        {
            if(light == "d" && ex < Sn_7Li)
                return true;
            else if(light == "p" && ex < Sn_8Li)
                return true;
            else
                return false;
        },
        {"Ex"});
    auto df_L1_unbound = df_L1.Filter(
        [&light, &Sn_7Li, &Sn_8Li](double ex)
        {
            if(light == "d" && ex > Sn_7Li)
                return true;
            else if(light == "p" && ex > Sn_8Li)
                return true;
            else
                return false;
        },
        {"Ex"});

    // Histograms for phisum distribution
    auto hPhiSum = df.Histo1D(
        {"hPhiSum", "Phi sum distribution;#phi_{light} - #phi_{heavy} [deg];Counts", 100, -180, 180}, "phiSum");
    auto hPhiSum_sil = df_sil.Histo1D(
        {"hPhiSum_sil",
         "Phi sum distribution for events with light in silicons;#phi_{light} - #phi_{heavy} [deg];Counts", 100, -180,
         180},
        "phiSum");
    auto hPhiSum_L1 = df_L1.Histo1D(
        {"hPhiSum_L1", "Phi sum distribution for events with light in L1;#phi_{light} - #phi_{heavy} [deg];Counts", 100,
         -180, 180},
        "phiSum");

    // Histograms for phi distributions
    auto hPhiBeam = df.Histo1D({"hPhiBeam", "Phi of beam;#phi_{beam} [deg];Counts", 100, -180, 180}, "phiBeam");
    auto hPhiLight =
        df.Histo1D({"hPhiLight", "Phi of light;#phi_{light} [deg];Counts", 100, -180, 180}, "MergerData.fPhiLight");
    auto hPhiHeavy =
        df.Histo1D({"hPhiHeavy", "Phi of heavy;#phi_{heavy} [deg];Counts", 100, -180, 180}, "MergerData.fPhiHeavy");

    // Histogramas for Phi debugging
    auto hPhiLightLegacy = df.Histo1D(
        {"hPhiLightLegacy", "Phi of light (legacy);#phi_{light} [deg];Counts", 100, -180, 180}, "phiLightLegacy");
    auto hPhiHeavyLegacy = df.Histo1D(
        {"hPhiHeavyLegacy", "Phi of heavy (legacy);#phi_{heavy} [deg];Counts", 100, -180, 180}, "phiHeavyLegacy");
    auto hPhiLightNew =
        df.Histo1D({"hPhiLightNew", "Phi of light (new);#phi_{light} [deg];Counts", 100, -180, 180}, "phiLight");
    auto hPhiHeavyNew =
        df.Histo1D({"hPhiHeavyNew", "Phi of heavy (new);#phi_{heavy} [deg];Counts", 100, -180, 180}, "phiHeavy");
    auto hPhiSumLegacy = df.Histo1D(
        {"hPhiSumLegacy", "Phi sum distribution (legacy);#phi_{light} - #phi_{heavy} [deg];Counts", 100, -280, 180},
        "phiSumLegacy");
    auto hPhiSumNew = df.Histo1D(
        {"hPhiSumNew", "Phi sum distribution (new);#phi_{light} - #phi_{heavy} [deg];Counts", 100, -280, 180},
        "phiSumNew");

    // Histograms for phi dist and phisum for bound and unbound events
    auto hPhiSum_sil_bound = df_sil_bound.Histo1D(
        {"hPhiSum_sil_bound",
         "Phi sum distribution for events with light in silicons and bound;#phi_{light} - #phi_{heavy} [deg];Counts",
         100, -180, 180},
        "phiSum");
    auto hPhiSum_sil_unbound = df_sil_unbound.Histo1D(
        {"hPhiSum_sil_unbound",
         "Phi sum distribution for events with light in silicons and unbound;#phi_{light} - #phi_{heavy} [deg];Counts",
         100, -180, 180},
        "phiSum");
    auto hPhiSum_L1_bound = df_L1_bound.Histo1D(
        {"hPhiSum_L1_bound",
         "Phi sum distribution for events with light in L1 and bound;#phi_{light} - #phi_{heavy} [deg];Counts", 100,
         -180, 180},
        "phiSum");
    auto hPhiSum_L1_unbound = df_L1_unbound.Histo1D(
        {"hPhiSum_L1_unbound",
         "Phi sum distribution for events with light in L1 and unbound;#phi_{light} - #phi_{heavy} [deg];Counts", 100,
         -180, 180},
        "phiSum");

    auto hPhiLight_sil_bound = df_sil_bound.Histo1D(
        {"hPhiLight_sil_bound", "Phi of light in silicons and bound;#phi_{light} [deg];Counts", 100, -180, 180},
        "MergerData.fPhiLight");
    auto hPhiLight_sil_unbound = df_sil_unbound.Histo1D(
        {"hPhiLight_sil_unbound", "Phi of light in silicons and unbound;#phi_{light} [deg];Counts", 100, -180, 180},
        "MergerData.fPhiLight");
    auto hPhiLight_L1_bound = df_L1_bound.Histo1D(
        {"hPhiLight_L1_bound", "Phi of light in L1 and bound;#phi_{light} [deg];Counts", 100, -180, 180},
        "MergerData.fPhiLight");
    auto hPhiLight_L1_unbound = df_L1_unbound.Histo1D(
        {"hPhiLight_L1_unbound", "Phi of light in L1 and unbound;#phi_{light} [deg];Counts", 100, -180, 180},
        "MergerData.fPhiLight");

    auto hPhiHeavy_sil_bound = df_sil_bound.Histo1D(
        {"hPhiHeavy_sil_bound", "Phi of heavy in silicons and bound;#phi_{heavy} [deg];Counts", 100, -180, 180},
        "MergerData.fPhiHeavy");
    auto hPhiHeavy_sil_unbound = df_sil_unbound.Histo1D(
        {"hPhiHeavy_sil_unbound", "Phi of heavy in silicons and unbound;#phi_{heavy} [deg];Counts", 100, -180, 180},
        "MergerData.fPhiHeavy");
    auto hPhiHeavy_L1_bound = df_L1_bound.Histo1D(
        {"hPhiHeavy_L1_bound", "Phi of heavy in L1 and bound;#phi_{heavy} [deg];Counts", 100, -180, 180},
        "MergerData.fPhiHeavy");
    auto hPhiHeavy_L1_unbound = df_L1_unbound.Histo1D(
        {"hPhiHeavy_L1_unbound", "Phi of heavy in L1 and unbound;#phi_{heavy} [deg];Counts", 100, -180, 180},
        "MergerData.fPhiHeavy");

    // Histograms for phisSum correlations
    // Phisum - Qtot correlation
    auto hPhiSum_Qtot_sil =
        df_sil.Histo2D({"hPhiSum_Qtot_sil", "Phi sum vs Qtot;#phi_{light} - #phi_{heavy} [deg];Q_{tot} [pC]", 100, -180,
                    180, 2000, 0, 3e5},
                   "phiSum", "MergerData.fLight.fQtotal");
    auto hPhiSum_Qtot_L1 =
        df_L1.Histo2D({"hPhiSum_Qtot_L1", "Phi sum vs Qtot;#phi_{light} - #phi_{heavy} [deg];Q_{tot} [pC]", 100, -180, 180,
                    2000, 0, 3e5},
                   "phiSum", "MergerData.fLight.fQtotal");

    // Phisum - Ex correlation
    auto hPhiSum_Ex_sil = df_sil.Histo2D(
        {"hPhiSum_Ex_sil", "Phi sum vs Ex;#phi_{light} - #phi_{heavy} [deg];E_{x} [MeV]", 100, -180, 180, 100, -5, 10},
        "phiSum", "Ex");
    auto hPhiSum_Ex_L1 = df_L1.Histo2D(
        {"hPhiSum_Ex_L1", "Phi sum vs Ex;#phi_{light} - #phi_{heavy} [deg];E_{x} [MeV]", 100, -180, 180, 100, -5, 10},
        "phiSum", "Ex");


    auto* cPhiSum = new TCanvas("cPhiSum", "Phi sum distribution", 800, 600);
    cPhiSum->Divide(3, 1);
    cPhiSum->cd(1);
    hPhiSum->DrawClone();
    cPhiSum->cd(2);
    hPhiSum_sil->DrawClone();
    cPhiSum->cd(3);
    hPhiSum_L1->DrawClone();

    auto* cPhi = new TCanvas("cPhi", "Phi distributions", 800, 600);
    cPhi->Divide(3, 1);
    cPhi->cd(1);
    hPhiBeam->DrawClone();
    cPhi->cd(2);
    hPhiLight->DrawClone();
    cPhi->cd(3);
    hPhiHeavy->DrawClone();

    auto* cPhiDebug = new TCanvas("cPhiDebug", "Phi distributions debug", 800, 600);
    cPhiDebug->Divide(2, 3);
    cPhiDebug->cd(1);
    hPhiLightLegacy->DrawClone();
    cPhiDebug->cd(2);
    hPhiLightNew->DrawClone();
    cPhiDebug->cd(3);
    hPhiHeavyLegacy->DrawClone();
    cPhiDebug->cd(4);
    hPhiHeavyNew->DrawClone();
    cPhiDebug->cd(5);
    hPhiSumLegacy->DrawClone();
    cPhiDebug->cd(6);
    hPhiSumNew->DrawClone();

    auto* cCorrelations = new TCanvas("cCorrelations", "Correlations", 800, 600);
    cCorrelations->Divide(2, 2);
    cCorrelations->cd(1);
    hPhiSum_Ex_sil->DrawClone("colz");
    cCorrelations->cd(2);
    hPhiSum_Qtot_sil->DrawClone("colz");
    cCorrelations->cd(3);
    hPhiSum_Ex_L1->DrawClone("colz");
    cCorrelations->cd(4);
    hPhiSum_Qtot_L1->DrawClone("colz");

    auto* cBoundUnbound = new TCanvas("cBoundUnbound", "Bound and unbound distributions", 800, 600);
    cBoundUnbound->Divide(2, 2);
    const std::vector<std::string> labels = {"Silicon - bound", "Silicon - unbound", "L1 - bound", "L1 - unbound"};
    auto drawGroup = [&](TVirtualPad* pad, std::vector<TH1*> hs, bool drawLegend = true)
    {
        const int colors[] = {kRed, kBlue, kGreen + 2, kOrange + 1};

        // 1) máximo de todos los histogramas del pad
        double ymax = 0;
        for(auto* h : hs)
            ymax = std::max(ymax, h->GetMaximum());

        // 2) dibujar con color y máximo ya fijados, y leyenda con los clones
        pad->cd();
        auto* leg = new TLegend(0.55, 0.65, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);

        for(size_t i = 0; i < hs.size(); ++i)
        {
            hs[i]->SetLineColor(colors[i]);
            hs[i]->SetMaximum(1.2 * ymax);
            auto* clone = (TH1*)hs[i]->DrawClone(i == 0 ? "" : "same");
            leg->AddEntry(clone, labels[i].c_str(), "l");
        }
        if(drawLegend)
            leg->Draw();
    };
    drawGroup(cBoundUnbound->cd(1), {hPhiSum_sil_bound.GetPtr(), hPhiSum_sil_unbound.GetPtr(),
                                     hPhiSum_L1_bound.GetPtr(), hPhiSum_L1_unbound.GetPtr()});
    drawGroup(cBoundUnbound->cd(2), {hPhiLight_sil_bound.GetPtr(), hPhiLight_sil_unbound.GetPtr(),
                                     hPhiLight_L1_bound.GetPtr(), hPhiLight_L1_unbound.GetPtr()});
    drawGroup(cBoundUnbound->cd(3), {hPhiHeavy_sil_bound.GetPtr(), hPhiHeavy_sil_unbound.GetPtr(),
                                     hPhiHeavy_L1_bound.GetPtr(), hPhiHeavy_L1_unbound.GetPtr()});
}