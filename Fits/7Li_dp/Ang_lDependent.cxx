#include "ActMergerData.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TColor.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"

#include "AngComparator.h"
#include "AngDifferentialXS.h"
#include "AngFitter.h"
#include "AngGlobals.h"
#include "AngIntervals.h"
#include "FitInterface.h"
#include "Interpolators.h"
#include "PhysExperiment.h"

#include <string>
#include <vector>

#include "../../PrettyStyle.C"
#include "../Histos.h"

void Ang_lDependent(bool isLab = false)
{
    PrettyStyle(false);
    if(isLab)
        Angular::ToggleIsLab();

    ROOT::EnableImplicitMT();

    ROOT::RDataFrame df {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_F_7Li_d_p_filtered.root"};
    auto def {df.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == true; },
                        {"MergerData"})}; // only silicons, == false is for L1 events
    ROOT::RDataFrame phase {"SimulationTTree", "../../Simulation/Outputs/7Li/2H_1H_TRIUMF_Eex_0.000_nPS_1_pPS_0.root"};
    // Book histograms
    auto hEx {def.Histo1D(S2384Fit::Exdp_7Li, "Ex")};
    auto hCM {def.Histo2D({"hCM", "CM;#theta_{CM};E [MeV]", 300, 0, 120, 300, 0, 60}, "ThetaCM", "EVertex")};

    // Init intervals
    double thetaMin {32};
    double thetaMax {80};
    double thetaStep {5};
    Angular::Intervals ivs {thetaMin, thetaMax, S2384Fit::Exdp_7Li, thetaStep, 1};
    def.Foreach([&](double thetacm, double ex) { ivs.Fill(thetacm, ex); }, {"ThetaCM", "Ex"});
    phase.Foreach([&](double thetacm, double ex, double weight) { ivs.FillPS(0, thetacm, ex, weight); },
                  {"theta3CM", "Eex", "weight"});
    ivs.TreatPS(0, 0.2);
    ivs.Draw();

    // Init fitter
    Angular::Fitter fitter {&ivs};
    fitter.SetAllowFreeMean(false);
    // fitter.SetFreeMeanRange(0.1);
    fitter.Configure("./Outputs/fit_lDependent.root");
    double R {(std::pow(7, 1. / 3.) + std::pow(1, 1. / 3.)) * 1.25}; // interaction radius in fm, with r0 = 1.25 fm
    double mu {7 * 1. / (7 + 1) * 931.5};                            // reduced mass in MeV/c^2
    fitter.ApplyLambdaToModels([mu, R](Fitters::Model& m) { m.AddBWL(0, 2, 2.03262, mu, R); });
    fitter.Run();
    fitter.Draw();
    fitter.DrawCounts();

    // Interface
    Fitters::Interface inter;
    inter.Read("./Outputs/interface.root");
    auto peaks {inter.GetKeys()};

    // Efficiency
    Interpolators::Efficiency eff;
    eff.Add("g0", "../../Simulation/Outputs/7Li/2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0.root", isLab ? "effLabside" : "effCMside");
    eff.Add("g1", "../../Simulation/Outputs/7Li/2H_1H_TRIUMF_Eex_0.981_nPS_0_pPS_0.root", isLab ? "effLabside" : "effCMside");
    eff.Add("g2", "../../Simulation/Outputs/7Li/2H_1H_TRIUMF_Eex_2.255_nPS_0_pPS_0.root", isLab ? "effLabside" : "effCMside");
    eff.Add("v0", "../../Simulation/Outputs/7Li/2H_1H_TRIUMF_Eex_3.210_nPS_0_pPS_0.root", isLab ? "effLabside" : "effCMside");
    eff.Add("v1", "../../Simulation/Outputs/7Li/2H_1H_TRIUMF_Eex_5.400_nPS_0_pPS_0.root", isLab ? "effLabside" : "effCMside");
    eff.Add("v2", "../../Simulation/Outputs/7Li/2H_1H_TRIUMF_Eex_6.100_nPS_0_pPS_0.root", isLab ? "effLabside" : "effCMside");
    eff.Add("v3", "../../Simulation/Outputs/7Li/2H_1H_TRIUMF_Eex_6.530_nPS_0_pPS_0.root", isLab ? "effLabside" : "effCMside");
    eff.Add("v4", "../../Simulation/Outputs/7Li/2H_1H_TRIUMF_Eex_7.100_nPS_0_pPS_0.root", isLab ? "effLabside" : "effCMside");
    eff.Add("v5", "../../Simulation/Outputs/7Li/2H_1H_TRIUMF_Eex_7.100_nPS_0_pPS_0.root", isLab ? "effLabside" : "effCMside");
    // eff.Add("ps0", "../../Simulation/Outputs/7Li/2H_1H_TRIUMF_Eex_0.000_nPS_1_pPS_0.root", isLab ? "effLab" :
    // "effCM");
    //  eff.Add("g1", "./Inputs/effs/g1_7Li_dp_sil.root", "effCM");
    //   Draw to check is fine
    eff.Draw();

    // Set experiment info
    PhysUtils::Experiment exp {"../norm/7Li_norm.dat"};
    // And compute differential xs!
    Angular::DifferentialXS xs {&ivs, &fitter, &eff, &exp};
    xs.DoFor(peaks);
    xs.Write("./Outputs/", "xs_lDependent");

    // Plot
    // Comparators!
    for(const auto& peak : peaks)
        inter.AddAngularDistribution(peak, xs.Get(peak));
    inter.ReadCompConfig("./comps.conf");
    inter.FillComp();
    inter.FitComp();

    Angular::Comparator comp {"g.s", xs.Get("g0")};
    comp.Add("ADWA", "./Inputs/gs/21.gs");
    comp.Add("DA1p-Delaroche", "./Inputs/gs_DA1p_Delaroche/21.g0");
    comp.Add("Daehnik-Delaroche", "./Inputs/gs_Daehnik_Delaroche/21.g0");
    comp.Add("DA1pcorr-Delaroche", "./Inputs/gs_DA1pcorr_Delaroche/21.g0");
    Angular::Comparator comp1 {"1st Ex", xs.Get("g1")};
    comp1.Add("Daehnik-Delaroche 1st Ex", "./Inputs/g1_Daehnik_Delaroche/21.g1");
    comp1.Add("DA1pcorr-Delaroche 1st Ex", "./Inputs/g1_DA1pcorr_Delaroche/21.g1");
    Angular::Comparator comp2 {"2nd Ex", xs.Get("g2")};
    comp2.Add("Daehnik-Delaroche 2nd Ex", "./Inputs/g2_Daehnik_Delaroche/21.g2");
    comp2.Add("DA1pcorr-Delaroche 2nd Ex", "./Inputs/g2_DA1pcorr_Delaroche/21.g2");
    Angular::Comparator comp3 {"3rd Ex", xs.Get("v0")};
    comp3.Add("Daehnik-Delaroche 3rd Ex", "./Inputs/g3_Daehnik_Delaroche/21.g3");
    comp3.Add("DA1pcorr-Delaroche 3rd Ex", "./Inputs/g3_DA1pcorr_Delaroche/21.g3");
    comp.Fit();
    comp.Draw("gs", true);
    comp.DrawTheo();
    comp1.Fit();
    comp1.Draw("1st", true);
    comp2.Fit();
    comp2.Draw("2nd", true);
    comp2.DrawTheo();
    comp3.Fit();
    comp3.Draw("3rd", true);
    comp3.DrawTheo();

    auto* c0 {new TCanvas {"c0", "(d,p) canvas"}};
    c0->DivideSquare(2);
    c0->cd(1);
    hEx->DrawClone();
    c0->cd(2);
    hCM->DrawClone("colz");
}
