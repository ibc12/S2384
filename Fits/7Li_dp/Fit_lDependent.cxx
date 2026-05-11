#include "ActMergerData.h"

#include "ROOT/RDataFrame.hxx"

#include "TROOT.h"

#include "FitInterface.h"
#include "FitModel.h"
#include "FitUtils.h"
#include "Interpolators.h"

#include <cmath>
#include <string>
#include <vector>

#include "../../PrettyStyle.C"
#include "../Histos.h"
void Fit_lDependent()
{
    PrettyStyle();
    ROOT::EnableImplicitMT();

    // Analysis
    ROOT::RDataFrame df {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_F_7Li_d_p_filtered.root"};
    auto def {df.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == true; }, {"MergerData"})};
    // Ex
    auto hEx {def.Histo1D(S2384Fit::Exdp_7Li, "Ex")};
    // Phase space
    ROOT::RDataFrame phase {"SimulationTTree", "../../Simulation/Outputs/7Li/2H_1H_TRIUMF_Eex_0.000_nPS_1_pPS_0.root"};
    auto hPS {phase.Histo1D(S2384Fit::Exdp_7Li, "Eex", "weight")};

    // Sigmas
    Interpolators::Sigmas sigmas;
    sigmas.Read("../../Simulation/Outputs/7Li/sigmas_7Li_2H_1H.root");

    // Interface to fit
    Fitters::Interface inter;
    double sigma {0.14}; // common init sigma for all
    double gamma {0.05}; // common init gamma for all voigts
    inter.AddState("g0", {100, 0, sigma});
    inter.AddState("g1", {30, 0.98, sigma});
    inter.AddState("g2", {60, 2.2, sigma});
    inter.AddState("v0", {20, 3.2, sigma, 1});
    // inter.AddState("v1", {15, 5.4, sigma, 0.65});
    // inter.AddState("v2", {15, 6.1, sigma, 1});
    // inter.AddState("v3", {10, 6.5, sigma, 0.035});
    // inter.AddState("v4", {10, 7.1, sigma, 0.4});
    inter.AddState("ps0", {1e-6}, "ps0");
    inter.EndAddingStates();
    inter.EvalSigma(sigmas.GetGraph());
    inter.SetFixAll(2, true); // fix all sigmas
    // inter.SetFixAll(3, true); // fix all gammas
    // inter.SetBoundsAll(2, {0.05, 0.3}); // sigma bounds
    inter.SetBounds("g0", 0, {10, 200});
    inter.SetBounds("g1", 0, {10, 100});
    inter.SetBounds("g2", 0, {10, 100});
    inter.SetBounds("v0", 3, {0.5, 1.5}); // gamma bounds
    // inter.SetBounds("v1", 3, {0.3, 1});
    // inter.SetBounds("v2", 3, {0.5, 1.5});
    // inter.SetBounds("v3", 3, {0.01, 0.1});
    // inter.SetBounds("v4", 3, {0.2, 0.8});

    // Save to be used later
    inter.Write("./Outputs/interface.root");

    // Model
    Fitters::Model model {inter.GetNGauss(), inter.GetNVoigt(), {*hPS}};
    double R {(std::pow(7, 1. / 3.) + std::pow(1, 1. / 3.)) * 1.25}; // interaction radius in fm, with r0 = 1.25 fm
    double mu {7 * 1. / (7 + 1) * 931.5}; // reduced mass in MeV/c^2
    model.AddGammaL(0, 2, 2.03262, mu, R);

    // Fitting range
    double exmin {-2};
    double exmax {9};

    // Run!
    Fitters::RunFit(hEx.GetPtr(), exmin, exmax, model, inter.GetInitial(), inter.GetBounds(), inter.GetFixed(),
                    ("./Outputs/fit_lDependent.root"), "7Li(d,p) fit",
                    {{"g0", "g.s"}, {"g1", "1st ex"}, {"g2", "2nd ex"}, {"ps0", "1-n phase"}}, false);
}
