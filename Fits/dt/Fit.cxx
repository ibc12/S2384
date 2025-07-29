#include "ROOT/RDataFrame.hxx"

#include "TROOT.h"

#include "FitInterface.h"
#include "FitModel.h"
#include "FitUtils.h"

#include <string>
#include <vector>

#include "../Histos.h"
void Fit()
{
    ROOT::EnableImplicitMT();

    // Analysis
    ROOT::RDataFrame df {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_d_t.root"};
    // Ex
    auto hEx {df.Histo1D(S2384Fit::Exdt, "Ex")};

    // Interface to fit
    Fitters::Interface inter;
    double sigma {0.3}; // common init sigma for all
    inter.AddState("g0", {400, 0, sigma}, "gs");
    inter.EndAddingStates();
    inter.SetBounds("g0", 1, {-0.2, 0.4});
    inter.SetBounds("g0", 2, {0.1, 4});
    // Save to be used later
    inter.Write("./Outputs/interface.root");

    // Model
    Fitters::Model model {inter.GetNGauss(), inter.GetNVoigt(), {}};

    // Fitting range
    double exmin {-2};
    double exmax {1.5};

    // Run!
    Fitters::RunFit(hEx.GetPtr(), exmin, exmax, model, inter.GetInitial(), inter.GetBounds(), inter.GetFixed(),
                    ("./Outputs/fit.root"), "11Li(d,t) fit", {{"g0", "g.s"}}, false);
}
