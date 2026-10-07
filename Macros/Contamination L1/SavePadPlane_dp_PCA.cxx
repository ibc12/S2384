#include "ActMergerData.h"
#include "ActParticle.h"
#include "ActSRIM.h"

#include "ROOT/RDataFrame.hxx"

#include "TMath.h"
#include "TString.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// Transverse momentum (w.r.t. beam axis) of the light particle. Returns -999 if unphysical
double LightTransverse(double thL, double ELightKin, double mLight)
{
    double ELight {ELightKin + mLight};
    double pL2 {ELight * ELight - mLight * mLight};
    if(!(pL2 >= 0)) // also catches NaN
        return -999.;
    return std::sqrt(pL2) * std::sin(thL);
}

// Transverse momentum of the heavy particle. Its momentum comes from energy conservation
// (heavy mass + Ex) and the transverse component uses the MEASURED heavy angle. Returns -999 if unphysical
double HeavyTransverse(double thH, double EBeam, double ELightKin, double Ex, double mBeam, double mTarget,
                       double mLight, double mHeavy)
{
    double EHeavy {(EBeam + mBeam) + mTarget - (ELightKin + mLight)};
    double mEx {mHeavy + Ex};
    double pH2 {EHeavy * EHeavy - mEx * mEx};
    if(!(pH2 >= 0))
        return -999.;
    return std::sqrt(pH2) * std::sin(thH);
}

// Saves, only for events whose light particle stops in the pad plane (no silicon energy),
// a flat tree with PID/angle variables, Ex and the transverse momentum of light and heavy
// under both the (d,p) and the (d,d) hypotheses.
void SavePadPlane_dp_PCA()
{
    std::string beam {"7Li"};
    std::string light {"p"};

    ROOT::RDataFrame dfini {"Final_Tree", TString::Format("../../PostAnalysis/Outputs/tree_ex_%s_d_%s_filtered.root",
                                                          beam.c_str(), light.c_str())};
    TString outName {TString::Format("./Outputs/tree_padplane_%s_d_%s.root", beam.c_str(), light.c_str())};

    // Masses [MeV]
    const double mBeam {ActPhysics::Particle {"7Li"}.GetMass()};
    const double mD {ActPhysics::Particle {"d"}.GetMass()};
    const double mP {ActPhysics::Particle {"p"}.GetMass()};
    const double m8Li {ActPhysics::Particle {"8Li"}.GetMass()};
    const double m7Li {mBeam};

    // SRIM for the deuteron
    auto* srim {new ActPhysics::SRIM};
    // NOTE: path relative to this macro. Adjust if needed
    srim->ReadTable("d", "../../Calibrations/SRIM/2H_900mb_CF4_95-5.txt");

    auto dfPad =
        dfini
            // Only particles stopped in the pad plane (L1 trigger branch of your EVertex calculation)
            .Filter([](const ActRoot::MergerData& d) { return !d.fLight.IsFilled(); }, {"MergerData"}, "Pad plane")
            // ---- Variables to save (new names, to avoid clashing with columns already in the tree)
            // NOTE: fQave is read as d.fQave. If it lives inside fLight, change it to d.fLight.fQave
            .Define("Qave", [](const ActRoot::MergerData& d) { return static_cast<double>(d.fQave); }, {"MergerData"})
            .Define("Qtotal", [](const ActRoot::MergerData& d) { return static_cast<double>(d.fLight.fQtotal); },
                    {"MergerData"})
            .Define("TL", [](const ActRoot::MergerData& d) { return static_cast<double>(d.fLight.fTL); },
                    {"MergerData"})
            .Define("ThetaLight", [](const ActRoot::MergerData& d) { return static_cast<double>(d.fThetaLight); },
                    {"MergerData"})
            .Define("ThetaHeavy", [](const ActRoot::MergerData& d) { return static_cast<double>(d.fThetaHeavy); },
                    {"MergerData"})
            .Define("PhiLight", [](const ActRoot::MergerData& d) { return static_cast<double>(d.fPhiLight); },
                    {"MergerData"})
            .Define("PhiHeavy", [](const ActRoot::MergerData& d) { return static_cast<double>(d.fPhiHeavy); },
                    {"MergerData"})
            // ---- Auxiliary (not saved): angles in rad and light energy as a deuteron
            .Define("thL_rad", [](double th) { return th * TMath::DegToRad(); }, {"ThetaLight"})
            .Define("thH_rad", [](double th) { return th * TMath::DegToRad(); }, {"ThetaHeavy"})
            .Define("EVertex_d", [srim](const ActRoot::MergerData& d)
                    { return srim->EvalEnergy("d", d.fLight.fTL); }, // stopped in pad plane: range only
                    {"MergerData"})
            // ---- Transverse momentum of each particle, under each hypothesis
            // (d,p): light = p (EVertex from the tree), heavy = 8Li + Ex (from the tree)
            .Define("pTLight_dp", [=](double thL, double EVertex) { return LightTransverse(thL, EVertex, mP); },
                    {"thL_rad", "EVertex"})
            .Define("pTHeavy_dp", [=](double thH, double EBeam, double EVertex, double Ex)
                    { return HeavyTransverse(thH, EBeam, EVertex, Ex, mBeam, mD, mP, m8Li); },
                    {"thH_rad", "EBeam", "EVertex", "Ex"})
            // (d,d): light = d (EVertex_d from deuteron SRIM), heavy = 7Li g.s. (Ex = 0)
            .Define("pTLight_dd", [=](double thL, double EVertex_d) { return LightTransverse(thL, EVertex_d, mD); },
                    {"thL_rad", "EVertex_d"})
            .Define("pTHeavy_dd", [=](double thH, double EBeam, double EVertex_d)
                    { return HeavyTransverse(thH, EBeam, EVertex_d, 0., mBeam, mD, mD, m7Li); },
                    {"thH_rad", "EBeam", "EVertex_d"})
            .Define("phiSum",
                    [](ActRoot::MergerData& m)
                    {
                        double substraction {0};
                        if(m.fPhiLight > m.fPhiHeavy)
                            substraction += (m.fPhiLight - m.fPhiHeavy);
                        else
                            substraction += (m.fPhiHeavy - m.fPhiLight);
                        return substraction - 180;
                    },
                    {"MergerData"});

    std::vector<std::string> cols {"Qave", "Qtotal",     "ThetaLight", "TL",         "phiSum",    "ThetaHeavy",
                                   "Ex",   "pTLight_dp", "pTHeavy_dp", "pTLight_dd", "pTHeavy_dd"};

    // Booked before the Snapshot, so it is filled in the same event loop
    auto nEvents {dfPad.Count()};
    dfPad.Filter([](double ex) { return ex > 2; }, {"Ex"}, "Pad plane").Snapshot("PadPlane_Tree", outName.Data(), cols);

    std::cout << "Saved " << *nEvents << " pad plane events in " << outName << '\n';
}