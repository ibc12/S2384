#include "ROOT/RDataFrame.hxx"

#include "TROOT.h"

#include "FitInterface.h"
#include "FitModel.h"
#include "FitUtils.h"

#include "ActMergerData.h"

#include <string>
#include <vector>

#include "../Histos.h"
void Fit()
{
    ROOT::EnableImplicitMT();

    // Analysis
    ROOT::RDataFrame df {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_7Li_d_p_filtered.root"};
    auto def {df.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == true; }, {"MergerData"})};
    // Ex
    auto hEx {def.Histo1D(S2384Fit::Exdp_7Li, "Ex")};
    // Phase space 
    // ROOT::RDataFrame phase {"SimulationTTree", "../../Simulation/Outputs/7Li/2H_1H_TRIUMF_Eex_0.000_nPS_1_pPS_0.root"};
    // auto hPS {phase.Histo1D(S2384Fit::Exdp_7Li, "Eex", "weight")};

    // Interface to fit
    Fitters::Interface inter;
    double sigma {0.14}; // common init sigma for all
    double gamma {0.05}; // common init gamma for all voigts
    inter.AddState("g0", {100, 0, sigma});
    // inter.AddState("g1", {30, 0.9, sigma});
    // inter.AddState("g2", {60, 2.0, sigma});
    // inter.AddState("v0", {20, 3.2, sigma, gamma});
    // inter.AddState("v1", {15, 5.4, sigma, gamma});
    // inter.AddState("v2", {10, 6.5, sigma, gamma});
    // inter.AddState("v3", {10, 7.1, sigma, gamma});
    // inter.AddState("ps0", {1.}, "ps0");
    inter.EndAddingStates();
    inter.SetFixAll(2, true); // fix all sigmas
    //inter.SetFixAll(3, true); // fix all gammas
    // Save to be used later
    inter.Write("./Outputs/interface.root");

    // Model
    Fitters::Model model {inter.GetNGauss(), inter.GetNVoigt()};

    // Fitting range
    double exmin {-2};
    double exmax {8.5};

    // Run!
    Fitters::RunFit(hEx.GetPtr(), exmin, exmax, model, inter.GetInitial(), inter.GetBounds(), inter.GetFixed(),
                    ("./Outputs/fit.root"), "7Li(d,p) fit", {{"g0", "g.s"}, {"g1", "1st ex"}, {"ps0", "1-n phase"}}, false);
}
