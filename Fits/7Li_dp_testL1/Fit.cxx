#include "ActMergerData.h"

#include "ROOT/RDataFrame.hxx"

#include "TROOT.h"

#include "FitInterface.h"
#include "FitModel.h"
#include "FitUtils.h"
#include "Interpolators.h"

#include <string>
#include <vector>

#include "../../PrettyStyle.C"
#include "../Histos.h"
void Fit()
{
    PrettyStyle(false);
    ROOT::EnableImplicitMT();

    // Analysis
    ROOT::RDataFrame df {"Final_Tree", "./Inputs/tree_ex_F_7Li_d_p_filtered.root"};
    auto def {df.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; },
                        {"MergerData"})}; // == false is for L1 events
    // Ex
    auto hEx {def.Histo1D(S2384Fit::Exdp_7Li, "Ex")};
    // Phase space
    // ROOT::RDataFrame phase {"SimulationTTree",
    // "../../../Simulations/TRIUNF/11Li(d,p)/Outputs/7.5MeV/2H_1H_TRIUMF_Eex_0.000_nPS_1_pPS_0_silspecs_spacer_7Li.root"};
    // auto hPS {phase.Histo1D(S2384Fit::Exdp_7Li, "Eex", "weight")};

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
    double sigma_g0 {0.159188}; // given by simu
    double sigma_g1 {0.167294};
    inter.AddState("g0", {240, 0, sigma_g0});
    inter.AddState("g1", {30, 0.981, sigma_g1});
    inter.EndAddingStates();

    inter.EvalSigma(sigmas.GetGraph());
    // inter.SetFix("g1", 2, true); // fix all sigmas
    // inter.SetFix("g0", 2, true); // fix g.s. sigma
    //  inter.SetBounds("g1", 1, {0.4, 0.55});
    //  Save to be used later
    //  inter.SetFix("g1", 1, true);
    inter.Write("./Outputs/interface.root");

    // Model
    Fitters::Model model {inter.GetNGauss(), inter.GetNVoigt(), {}};

    // Fitting range
    double exmin {-1};
    double exmax {2};
    // double exmin {-0.8};
    // double exmax {1};

    // Run!
    Fitters::RunFit(hEx.GetPtr(), exmin, exmax, model, inter.GetInitial(), inter.GetBounds(), inter.GetFixed(),
                    ("./Outputs/fit.root"), "7Li(d,p) fit", {{"g0", "g.s"}, {"g1", "1st"}}, false);
}
