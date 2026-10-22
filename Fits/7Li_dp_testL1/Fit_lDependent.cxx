#include "ActMergerData.h"

#include "ROOT/RDataFrame.hxx"

#include "TROOT.h"

#include "FitInterface.h"
#include "FitModel.h"
#include "FitRunner.h"
#include "FitUtils.h"
#include "Interpolators.h"

#include <cmath>
#include <string>
#include <vector>

#include "../../PrettyStyle.C"
#include "../Histos.h"
void Fit_lDependent()
{
    PrettyStyle(false, true);
    ROOT::EnableImplicitMT();

    // Analysis
    ROOT::RDataFrame df {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_F_7Li_d_p_filtered.root"};
    auto def {df.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"})};
    // auto def {df.Filter([](ActRoot::MergerData& m) { return (m.fLight.IsFilled() && m.fHeavy.GetNLayers() == 2) ==
    // true; }, {"MergerData"})}; Ex
    auto hEx {def.Histo1D(S2384Fit::Exdp_7Li, "Ex")};
    // Phase space
    ROOT::RDataFrame phase {"SimulationTTree", "../../Simulation/Outputs/7Li/test_rpCut/"
                                               "2H_1H_TRIUMF_Eex_0.000_nPS_1_pPS_0_L1_240rpCut_paralelization.root"};
    auto hPS {phase.Histo1D(S2384Fit::Exdp_7Li, "Eex", "weight")};
    hPS->SetNameTitle("hPS", "1n PS");
    // Histogram of contamination
    ROOT::RDataFrame phase2 {"Final_Tree", "../../Macros/SimulationExp Comparison/Outputs/tree_ex_F_7Li_d_as_p.root"};
    auto hPS2 {phase2.Histo1D(S2384Fit::Exdp_7Li, "Ex_d_as_p")};
    hPS2->SetNameTitle("hPS2", "1n PS2");

    // Sigmas (how to for L1 ?¿)
    Interpolators::Sigmas sigmas;
    // sigmas.Read("../../Simulation/Outputs/7Li/test_ang_straggling/sigmas_7Li_2H_2H_1AngStr.root", "grSigma");
    // sigmas.Read("../../Simulation/Outputs/7Li/test_ang_straggling/sigmas_7Li_2H_2H_1-25AngStr.root", "grSigma");
    // sigmas.Read("../../Simulation/Outputs/7Li/test_ang_straggling/sigmas_7Li_2H_2H_1-3AngStr.root", "grsigmas");
    // sigmas.Read("../../Simulation/Outputs/7Li/test_ang_straggling/sigmas_7Li_2H_2H_1-35AngStr.root", "grsigmas");
    sigmas.Read("../../Simulation/Outputs/7Li/test_ang_straggling/sigmas_7Li_2H_2H_1-4AngStr.root", "grsigmas");
    // sigmas.Read("../../Simulation/Outputs/7Li/test_ang_straggling/sigmas_7Li_2H_2H_1-45AngStr.root", "grsigmas");
    // sigmas.Read("../../Simulation/Outputs/7Li/test_ang_straggling/sigmas_7Li_2H_2H_1-5AngStr.root", "grSigma");
    // sigmas.Read("../../Simulation/Outputs/7Li/test_ang_straggling/sigmas_7Li_2H_2H_1-65AngStr.root", "grSigma");
    // sigmas.Read("../../Simulation/Outputs/7Li/test_ang_straggling/sigmas_7Li_2H_2H_1-75AngStr.root", "grSigma");

    // Interface to fit
    Fitters::Interface inter;
    double sigma_g0 {0.1378}; // given by simu
    double sigma_g1 {0.167294};
    inter.AddState("g0", {240, 0, sigma_g0});
    inter.AddState("g1", {30, 0.981, sigma_g1});
    double sigma {0.1};
    inter.AddState("v0", {60, 2.2, sigma, 0.1});
    inter.AddState("v1", {20, 3.3, sigma, 1});
    inter.AddState("v2", {15, 4.1, sigma, 0.65});
    inter.AddState("v3", {15, 5.4, sigma, 0.65});
    inter.AddState("v4", {15, 6.1, sigma, 1});
    inter.AddState("v5", {10, 6.5, sigma, 0.035});
    inter.AddState("v6", {5, 7.1, sigma, 0.4});
    inter.AddState("v7", {5, 7.8, sigma, 0.4});
    inter.AddState("v8", {5, 8.6, sigma, 0.4});
    inter.AddState("v9", {5, 2.7, sigma, 0.4});
    inter.AddState("ps0", {0.005});
    inter.AddState("ps1", {0.005});
    inter.EndAddingStates();
    inter.EvalSigma(sigmas.GetGraph());
    inter.SetFixAll(2, true); // fix all sigmas
    inter.SetFixAll(2, true);
    inter.SetBounds("v0", 3, {0.005, 0.3}); // gamma bounds
    inter.SetBounds("v1", 3, {0.5, 1.5});   // gamma bounds
    inter.SetBounds("v2", 3, {0.01, 1});
    inter.SetBounds("v3", 3, {0.3, 1});
    inter.SetBounds("v4", 3, {0.5, 1.5});
    inter.SetBounds("v5", 3, {0.01, 0.1});
    inter.SetBounds("v6", 3, {0.1, 0.8});
    inter.SetBounds("v7", 3, {0.1, 0.8});
    inter.SetBounds("v8", 3, {0.1, 0.8});
    inter.SetBounds("v9", 3, {0.1, 0.8});
    inter.SetOffsetMeanBounds(0.3);
    // inter.SetFixAll(2, true);
    // inter.EvalSigma(sigmas.GetGraph());
    // inter.SetFix("g1", 2, true); // fix all sigmas
    // inter.SetFix("g0", 2, true); // fix g.s. sigma
    //  inter.SetBounds("g1", 1, {0.4, 0.55});
    //  Save to be used later
    //  inter.SetFix("g1", 1, true);
    inter.Write("./Outputs/interface_lDependent.root");

    // Model
    Fitters::Model model {inter.GetNGauss(), inter.GetNVoigt(), {*hPS, *hPS2}};

    // Fitting range
    double exmin {-1};
    double exmax {10};
    // double exmin {-0.8};
    // double exmax {1};

    double R {(std::pow(7, 1. / 3.) + std::pow(1, 1. / 3.)) * 1.25}; // interaction radius in fm, with r0 = 1.25 fm
    double mu {7 * 1. / (7 + 1) * 931.5};                            // reduced mass in MeV/c^2
    model.AddBWL(0, 1, 2.03262, mu, R);
    model.AddBWL(1, 1, 2.03262, mu, R);
    model.AddBWL(2, 1, 2.03262, mu, R);
    model.AddBWL(3, 2, 2.03262, mu, R);
    model.AddBWL(4, 2, 2.03262, mu, R);
    model.AddBWL(5, 1, 2.03262, mu, R);
    model.AddBWL(6, 2, 2.03262, mu, R);
    model.AddBWL(7, 2, 2.03262, mu, R);
    model.AddBWL(8, 2, 2.03262, mu, R);
    model.AddBWL(9, 1, 2.03262, mu, R);

    // Run!
    Fitters::RunFit(hEx.GetPtr(), exmin, exmax, model, inter.GetInitial(), inter.GetBounds(), inter.GetFixed(),
                    ("./Outputs/fit_lDependent.root"), "7Li(d,p) fit", {{"g0", "g.s"}, {"g1", "1st"}}, false);
}
