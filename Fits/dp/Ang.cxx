#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TROOT.h"
#include "TString.h"

#include "AngComparator.h"
#include "AngDifferentialXS.h"
#include "AngFitter.h"
#include "AngGlobals.h"
#include "AngIntervals.h"
#include "FitInterface.h"
#include "Interpolators.h"
#include "PhysExperiment.h"

#include "ActMergerData.h"

#include <string>
#include <vector>
#include "TStyle.h"
#include "TColor.h"

#include "../Histos.h"

void Ang()
{

    ROOT::EnableImplicitMT();

    ROOT::RDataFrame df {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_11Li_d_d.root"};
    auto def {df.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == true; }, {"MergerData"})}; // only silicons, == false is for L1 events

    // Book histograms
    auto hEx {def.Histo1D(S2384Fit::Exdd, "Ex")};
    auto hCM {def.Histo2D({"hCM", "CM;#theta_{CM};E [MeV]", 300, 0, 120, 300, 0, 60}, "ThetaCM", "EVertex")};

    // Init intervals
    double thetaMin {30};
    double thetaMax {60};
    double thetaStep {2};
    Angular::Intervals ivs {thetaMin, thetaMax, S2384Fit::Exdd, thetaStep, 0};
    def.Foreach([&](double thetacm, double ex) { ivs.Fill(thetacm, ex); }, {"ThetaCM", "Ex"});
    ivs.Draw();

    // Init fitter
    Angular::Fitter fitter {&ivs};
    fitter.SetAllowFreeMean(false);
    // fitter.SetFreeMeanRange(0.1);
    fitter.Configure("./Outputs/fit.root");
    fitter.Run();
    fitter.Draw();
    fitter.DrawCounts();

    // Interface
    Fitters::Interface inter;
    inter.Read("./Outputs/interface.root");
    auto peaks {inter.GetKeys()};

    // Efficiency
    Interpolators::Efficiency eff;
    eff.Add("g0", "./Inputs/effs/eff_11Li_dd_latsil.root", "effCM");
    // Draw to check is fine
    eff.Draw();

    // Set experiment info
    PhysUtils::Experiment exp {4.126e19 * 25.6, 9.87839e8, 1};
    // And compute differential xs!
    Angular::DifferentialXS xs {&ivs, &fitter, &eff, &exp};
    xs.DoFor(peaks);
    xs.Write("./Outputs/");

    // Plot
    Angular::Comparator comp {"g.s", xs.Get("g0")};
    comp.Add("Haixia", "./Inputs/gsH/fort.201");
    comp.Fit();
    comp.Add("Daehnick", "./Inputs/gsD/fort.201");
    comp.Fit();
    comp.Add("DA1p", "./Inputs/gsDA1p/fort.201");
    comp.Fit();
    comp.Draw("", true);

    auto* c0 {new TCanvas {"c0", "(d,d) canvas"}};
    c0->DivideSquare(2);
    c0->cd(1);
    hEx->DrawClone();
    c0->cd(2);
    hCM->DrawClone("colz");
}
