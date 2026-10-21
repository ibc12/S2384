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
    PrettyStyle(false, false);
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
    inter.AddState("ps0", {0.005});
    inter.AddState("ps1", {0.005});

    inter.EndAddingStates();
    // inter.SetFixAll(2, true);
    // inter.EvalSigma(sigmas.GetGraph());
    // inter.SetFix("g1", 2, true); // fix all sigmas
    // inter.SetFix("g0", 2, true); // fix g.s. sigma
    //  inter.SetBounds("g1", 1, {0.4, 0.55});
    //  Save to be used later
    //  inter.SetFix("g1", 1, true);
    inter.Write("./Outputs/interface.root");

    // Model
    Fitters::Model model {inter.GetNGauss(), inter.GetNVoigt(), {*hPS, *hPS2}};

    // Fitting range
    double exmin {-1};
    double exmax {10};
    // double exmin {-0.8};
    // double exmax {1};

    // Run!
    Fitters::RunFit(hEx.GetPtr(), exmin, exmax, model, inter.GetInitial(), inter.GetBounds(), inter.GetFixed(),
                    ("./Outputs/fit.root"), "7Li(d,p) fit", {{"g0", "g.s"}, {"g1", "1st"}}, false);
}
