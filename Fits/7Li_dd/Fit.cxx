#include "ROOT/RDataFrame.hxx"

#include "TROOT.h"

#include "FitInterface.h"
#include "FitModel.h"
#include "FitUtils.h"

#include <string>
#include <vector>

#include "ActMergerData.h"

#include "../Histos.h"
void Fit()
{
    ROOT::EnableImplicitMT();

    // Analysis
    ROOT::RDataFrame df {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_7Li_d_d.root"};
    auto def {df.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == true; }, {"MergerData"})}; // only silicons, == false is for L1 events
    // Ex
    auto hEx {def.Histo1D(S2384Fit::Exdd_7Li, "Ex")};
    // Phase space 
    // ROOT::RDataFrame phase {"SimulationTTree", "../../../Simulations/TRIUNF/11Li(d,p)/Outputs/7.5MeV/2H_1H_TRIUMF_Eex_0.000_nPS_1_pPS_0_silspecs_spacer_7Li.root"};
    // auto hPS {phase.Histo1D(S2384Fit::Exdp_7Li, "Eex", "weight")};

    // Interface to fit
    Fitters::Interface inter;
    double sigma {0.12}; // common init sigma for all
    inter.AddState("g0", {400, 0, sigma});
    inter.AddState("g1", {30, 0.435, sigma});
    inter.EndAddingStates();
    inter.SetFixAll(2, true); // fix all sigmas
    // Save to be used later
    inter.Write("./Outputs/interface.root");

    // Model
    Fitters::Model model {inter.GetNGauss(), inter.GetNVoigt(), {}};

    // Fitting range
    double exmin {-2};
    double exmax {4};

    // Run!
    Fitters::RunFit(hEx.GetPtr(), exmin, exmax, model, inter.GetInitial(), inter.GetBounds(), inter.GetFixed(),
                    ("./Outputs/fit.root"), "7Li(d,d) fit", {{"g0", "g.s"}}, false);
}
