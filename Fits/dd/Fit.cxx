#include "ROOT/RDataFrame.hxx"

#include "TROOT.h"

#include "FitInterface.h"
#include "FitModel.h"
#include "FitUtils.h"

#include <string>
#include <vector>

#include "./Histos.h"
void Fit()
{
    ROOT::EnableImplicitMT();

    // Analysis
    ROOT::RDataFrame df {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_d_d.root"};
    // Ex
    auto hEx {df.Histo1D(S2384Fit::Exdd, "Ex")};

    // Interface to fit
    Fitters::Interface inter;
    double sigma {0.3}; // common init sigma for all
    inter.AddState("g0", {400, 0, sigma}, "0+1");
    inter.EndAddingStates();
    // Save to be used later
    inter.Write("./Outputs/interface.root");

    // Model
    Fitters::Model model {inter.GetNGauss(), inter.GetNVoigt(), {}};

    // Fitting range
    double exmin {-2};
    double exmax {1};

    // Run!
    Fitters::RunFit(hEx.GetPtr(), exmin, exmax, model, inter.GetInitial(), inter.GetBounds(), inter.GetFixed(),
                    ("./Outputs/fit.root"), "20O(d,d) fit", {{"g0", "g.s"}}, false);
}
