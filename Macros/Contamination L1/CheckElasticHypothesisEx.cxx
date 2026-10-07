#include "ActMergerData.h"
#include "ActParticle.h"
#include "ActSRIM.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TMath.h"
#include "TStyle.h"
#include "TString.h"

#include <cmath>
#include <string>

// Transverse momentum imbalance pT_light - pT_heavy (w.r.t. beam axis) for a given hypothesis.
// The heavy momentum comes from energy conservation, using the heavy mass (+ Ex) and the MEASURED heavy angle.
//   ELightKin: light kinetic energy at the vertex [MeV] (must match the light mass hypothesis)
// Returns -999 for unphysical kinematics (outside the plot range, it is not a cut on the data).
double TransverseImbalance(double thL, double thH, double EBeam, double ELightKin, double Ex, double mBeam,
                           double mTarget, double mLight, double mHeavy)
{
    double ELight {ELightKin + mLight};
    double pL2 {ELight * ELight - mLight * mLight};
    double EHeavy {(EBeam + mBeam) + mTarget - ELight};
    double mEx {mHeavy + Ex};
    double pH2 {EHeavy * EHeavy - mEx * mEx};
    if(!(pL2 >= 0) || !(pH2 >= 0)) // also catches NaN
        return -999.;
    return std::sqrt(pL2) * std::sin(thL) - std::sqrt(pH2) * std::sin(thH);
}

// Correlations of the transverse momentum imbalance computed assuming
//   dp: light = p, heavy = 8Li (+Ex from the tree), EVertex from the tree (proton SRIM)
//   dd: light = d, heavy = 7Li g.s. (Ex = 0),       EVertex_d from the deuteron SRIM
// against each other and against Ex and Qtot. No cuts are applied.
void CheckElasticHypothesisEx()
{
    std::string beam {"7Li"};
    std::string light {"p"};

    ROOT::RDataFrame dfini {"Final_Tree", TString::Format("../../PostAnalysis/Outputs/tree_ex_%s_d_%s_filtered.root",
                                                          beam.c_str(), light.c_str())};

    // Masses [MeV]
    const double mBeam {ActPhysics::Particle {"7Li"}.GetMass()};
    const double mD {ActPhysics::Particle {"d"}.GetMass()};
    const double mP {ActPhysics::Particle {"p"}.GetMass()};
    const double m8Li {ActPhysics::Particle {"8Li"}.GetMass()};
    const double m7Li {mBeam};

    // ---- SRIM for the deuteron
    auto* srim {new ActPhysics::SRIM};
    // NOTE: path relative to this macro (your snippet uses the one relative to PostAnalysis/). Adjust if needed
    srim->ReadTable("d", "../../Calibrations/SRIM/2H_900mb_CF4_95-5.txt");

    auto df =
        dfini
            .Define("thL", [](const ActRoot::MergerData& d) { return static_cast<double>(d.fThetaLight) * TMath::DegToRad(); },
                    {"MergerData"})
            .Define("thH", [](const ActRoot::MergerData& d) { return static_cast<double>(d.fThetaHeavy) * TMath::DegToRad(); },
                    {"MergerData"})
            // PID variables (assumed member names: adjust fQtot if it is called differently)
            .Define("Qtot", [](const ActRoot::MergerData& d) { return static_cast<double>(d.fLight.fQtotal); },
                    {"MergerData"})
            .Define("TL", [](const ActRoot::MergerData& d) { return static_cast<double>(d.fLight.fTL); },
                    {"MergerData"})
            // Light energy at the vertex as if the particle were a deuteron
            .Define("EVertex_d",
                    [srim](const ActRoot::MergerData& d)
                    {
                        double ret {};
                        if(d.fLight.IsFilled())
                            ret = srim->EvalInitialEnergy("d", d.fLight.fEs.front(), d.fLight.fTL);
                        else // L1 trigger
                            ret = srim->EvalEnergy("d", d.fLight.fTL);
                        return ret;
                    },
                    {"MergerData"})
            // Transverse imbalance under each hypothesis
            .Define("pTDiff_dp",
                    [=](double thL, double thH, double EBeam, double EVertex, double Ex)
                    { return TransverseImbalance(thL, thH, EBeam, EVertex, Ex, mBeam, mD, mP, m8Li); },
                    {"thL", "thH", "EBeam", "EVertex", "Ex"})
            .Define("pTDiff_dd",
                    [=](double thL, double thH, double EBeam, double EVertex_d)
                    { return TransverseImbalance(thL, thH, EBeam, EVertex_d, 0., mBeam, mD, mD, m7Li); },
                    {"thL", "thH", "EBeam", "EVertex_d"});

    // Qtot and TL ranges from the data (single extra pass). Hardcode them if you prefer
    auto qMin {df.Min<double>("Qtot")};
    auto qMax {df.Max<double>("Qtot")};
    auto tlMin {df.Min<double>("TL")};
    auto tlMax {df.Max<double>("TL")};
    const double q0 {*qMin}, q1 {*qMax * 1.05};
    const double tl0 {*tlMin}, tl1 {*tlMax * 1.05};

    // ---- Histograms
    const int nPT {150};
    const double ptMin {-250}, ptMax {250};
    const int nEx {120};
    const double exMin {-5}, exMax {15};
    const int nQ {120};

    const char* tdp {"p_{T}^{L} - p_{T}^{H} (d,p) [MeV/c]"};
    const char* tdd {"p_{T}^{L} - p_{T}^{H} (d,d) [MeV/c]"};

    auto hDpDd = df.Histo2D({"hDpDd", TString::Format(";%s;%s", tdp, tdd), nPT, ptMin, ptMax, nPT, ptMin, ptMax},
                            "pTDiff_dp", "pTDiff_dd");
    auto hDpEx = df.Histo2D({"hDpEx", TString::Format(";E_{x} [MeV];%s", tdp), nEx, exMin, exMax, nPT, ptMin, ptMax},
                            "Ex", "pTDiff_dp");
    auto hDdEx = df.Histo2D({"hDdEx", TString::Format(";E_{x} [MeV];%s", tdd), nEx, exMin, exMax, nPT, ptMin, ptMax},
                            "Ex", "pTDiff_dd");
    auto hDpQ = df.Histo2D({"hDpQ", TString::Format(";Q_{tot} [arb. u.];%s", tdp), nQ, q0, q1, nPT, ptMin, ptMax},
                           "Qtot", "pTDiff_dp");
    auto hDdQ = df.Histo2D({"hDdQ", TString::Format(";Q_{tot} [arb. u.];%s", tdd), nQ, q0, q1, nPT, ptMin, ptMax},
                           "Qtot", "pTDiff_dd");
    auto hQTL = df.Histo2D({"hQTL", ";TL [mm];Q_{tot} [arb. u.]", nQ, tl0, tl1, nQ, q0, q1}, "TL", "Qtot");

    // ---- Draw
    gStyle->SetOptStat(0);
    auto* c = new TCanvas("cTransvHyp", "Transverse momentum: dp vs dd hypotheses", 1600, 1000);
    c->Divide(3, 2);

    c->cd(1);
    hDpDd->DrawClone("colz");
    c->cd(2);
    hDpEx->DrawClone("colz");
    c->cd(3);
    hDdEx->DrawClone("colz");
    c->cd(4);
    hQTL->DrawClone("colz");
    c->cd(5);
    hDpQ->DrawClone("colz");
    c->cd(6);
    hDdQ->DrawClone("colz");

    // for(int i {1}; i <= 6; i++)
    //     c->cd(i)->SetLogz();
}