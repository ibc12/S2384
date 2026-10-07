#include "ActInputParser.h"
#include "ActMergerData.h"
#include "ActParticle.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TMath.h"
#include "TString.h"

#include <cmath>
#include <string>

// Compares the transverse momentum (w.r.t. the beam axis) of the measured light particle
// with the one the heavy MUST carry if the reaction were a pure two-body channel.
// Bound events:   pT_light - pT_heavy ~ 0
// Unbound events: undetected neutron carries transverse momentum -> distribution shifts / widens
void CheckTransversalMomentumConservation()
{
    std::string beam {"7Li"};
    std::string target {"d"};
    std::string light {"p"};
    std::string heavy {"8Li"};

    // Get df from root file
    ROOT::RDataFrame dfini {"Final_Tree", TString::Format("../../PostAnalysis/Outputs/tree_ex_%s_d_%s_filtered.root",
                                                          beam.c_str(), light.c_str())};

    double Sn_8Li {2.03}; // MeV. Ex above this -> unbound

    // Masses [MeV]
    ActPhysics::Particle pBeam {beam};
    ActPhysics::Particle pTarget {target};
    ActPhysics::Particle pLight {light};
    ActPhysics::Particle pHeavy {heavy};
    const double mBeam {pBeam.GetMass()};
    const double mTarget {pTarget.GetMass()};
    const double mLight {pLight.GetMass()};
    const double mHeavy {pHeavy.GetMass()};

    auto df = dfini
                  // Coplanarity: ~0 for bound events
                  .Define("phiSum", [](ActRoot::MergerData& m) { return std::abs(m.fPhiLight - m.fPhiHeavy) - 180.; },
                          {"MergerData"})
                  // Light momentum from its kinetic energy at the vertex
                  .Define("pLight",
                          [=](double EVertex)
                          {
                              double E {EVertex + mLight};
                              return std::sqrt(E * E - mLight * mLight);
                          },
                          {"EVertex"})
                  // Heavy momentum expected from energy conservation, using the excited mass (mHeavy + Ex)
                  .Define("pHeavy",
                          [=](double EBeam, double EVertex, double Ex)
                          {
                              double EH {(EBeam + mBeam) + mTarget - (EVertex + mLight)};
                              double mEx {mHeavy + Ex};
                              double arg {EH * EH - mEx * mEx};
                              if(arg < 0)
                                  return -1.; // unphysical kinematics
                              return std::sqrt(arg);
                          },
                          {"EBeam", "EVertex", "Ex"})
                  // Transverse components w.r.t. beam axis
                  .Define("pTLight", [](ActRoot::MergerData& m, double pLight)
                          { return pLight * std::sin(m.fThetaLight * TMath::DegToRad()); }, {"MergerData", "pLight"})
                  .Define("pTHeavy", [](ActRoot::MergerData& m, double pHeavy)
                          { return pHeavy * std::sin(m.fThetaHeavy * TMath::DegToRad()); }, {"MergerData", "pHeavy"})
                  // Imbalance and normalized imbalance
                  .Define("pTDiff", [](double pTL, double pTH) { return pTL - pTH; }, {"pTLight", "pTHeavy"})
                  .Define("pTRel", [](double diff, double pTL) { return pTL != 0 ? diff / pTL : -999.; },
                          {"pTDiff", "pTLight"});

    // Keep only physical events
    auto dfPhys = df.Filter("pHeavy > 0", "Physical kinematics");
    auto dfBound = dfPhys.Filter([=](double Ex) { return Ex < Sn_8Li; }, {"Ex"}, "Bound");
    auto dfUnbound = dfPhys.Filter([=](double Ex) { return Ex >= Sn_8Li; }, {"Ex"}, "Unbound");

    // Histograms
    auto hDiffB = dfBound.Histo1D({"hDiffB", "Bound;p_{T}^{L} - p_{T}^{H} [MeV/c];Counts", 150, -200, 200}, "pTDiff");
    auto hDiffU =
        dfUnbound.Histo1D({"hDiffU", "Unbound;p_{T}^{L} - p_{T}^{H} [MeV/c];Counts", 150, -200, 200}, "pTDiff");
    auto hRelB = dfBound.Histo1D({"hRelB", "Bound;(p_{T}^{L} - p_{T}^{H}) / p_{T}^{L};Counts", 150, -2, 2}, "pTRel");
    auto hRelU =
        dfUnbound.Histo1D({"hRelU", "Unbound;(p_{T}^{L} - p_{T}^{H}) / p_{T}^{L};Counts", 150, -2, 2}, "pTRel");
    auto hVsEx = dfPhys.Histo2D({"hVsEx", ";E_{x} [MeV];p_{T}^{L} - p_{T}^{H} [MeV/c]", 120, -5, 15, 150, -200, 200},
                                "Ex", "pTDiff");
    auto hPhiSum = dfPhys.Histo2D(
        {"hPhiSum", ";E_{x} [MeV];#phi_{L} - #phi_{H} - 180 [#circ]", 120, -5, 15, 150, -60, 60}, "Ex", "phiSum");

    // Draw
    auto* c = new TCanvas("cTransversal", "Transverse momentum conservation", 1400, 900);
    c->Divide(3, 2);
    c->cd(1);
    hDiffB->SetLineColor(kBlue);
    hDiffU->SetLineColor(kRed);
    hDiffB->DrawNormalized("hist");
    hDiffU->DrawNormalized("hist same");
    c->cd(2);
    hRelB->SetLineColor(kBlue);
    hRelU->SetLineColor(kRed);
    hRelB->DrawNormalized("hist");
    hRelU->DrawNormalized("hist same");
    c->cd(3);
    hVsEx->DrawClone("colz");
    c->cd(4);
    hPhiSum->DrawClone("colz");
    c->cd(5);
    hDiffB->DrawClone("hist");
    c->cd(6);
    hDiffU->DrawClone("hist");
}