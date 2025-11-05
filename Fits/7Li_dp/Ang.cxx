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

#include <string>
#include <vector>
#include "TStyle.h"
#include "TColor.h"

#include "../Histos.h"

void Ang()
{

    ROOT::EnableImplicitMT();

    ROOT::RDataFrame df {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_7Li_d_p_filtered.root"};

    // Book histograms
    auto hEx {df.Histo1D(S2384Fit::Exdp_7Li, "Ex")};
    auto hCM {df.Histo2D({"hCM", "CM;#theta_{CM};E [MeV]", 300, 0, 120, 300, 0, 60}, "ThetaCM", "EVertex")};

    // Init intervals
    double thetaMin {32};
    double thetaMax {63};
    double thetaStep {5};
    Angular::Intervals ivs {thetaMin, thetaMax, S2384Fit::Exdd, thetaStep, 0};
    df.Foreach([&](double thetacm, double ex) { ivs.Fill(thetacm, ex); }, {"ThetaCM", "Ex"});
    ivs.Draw();

    // Init fitter
    Angular::Fitter fitter {&ivs};
    fitter.SetAllowFreeMean(true);
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
    eff.Add("g0", "./Inputs/effs/eff_7Li_dp_latsil.root", "effCM");
    // Draw to check is fine
    eff.Draw();

    // Set experiment info
    PhysUtils::Experiment exp {4.126e19*25.6, 1304030, 300};
    // And compute differential xs!
    Angular::DifferentialXS xs {&ivs, &fitter, &eff, &exp};
    xs.DoFor(peaks);
    xs.Write("./Outputs/");

    // Plot
    Angular::Comparator comp {"g.s", xs.Get("g0")};
    // comp.Add("Delaroche", "./Inputs/gs_DA1p_Delaroche/21.g0");
    comp.Add("Delaroche", "./Inputs/gs_Daehnik_Delaroche/21.g0");
    comp.Fit();
    comp.Draw("", true);

    auto* c0 {new TCanvas {"c0", "(d,p) canvas"}};
    c0->DivideSquare(2);
    c0->cd(1);
    hEx->DrawClone();
    c0->cd(2);
    hCM->DrawClone("colz");
}
