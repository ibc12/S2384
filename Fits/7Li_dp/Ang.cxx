#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TROOT.h"
#include "TString.h"

#include "ActMergerData.h"

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
    auto def {df.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == true; }, {"MergerData"})}; // only silicons, == false is for L1 events
    // Book histograms
    auto hEx {def.Histo1D(S2384Fit::Exdp_7Li, "Ex")};
    auto hCM {def.Histo2D({"hCM", "CM;#theta_{CM};E [MeV]", 300, 0, 120, 300, 0, 60}, "ThetaCM", "EVertex")};

    // Init intervals
    double thetaMin {32};
    double thetaMax {63};
    double thetaStep {5};
    Angular::Intervals ivs {thetaMin, thetaMax, S2384Fit::Exdd, thetaStep, 0};
    def.Foreach([&](double thetacm, double ex) { ivs.Fill(thetacm, ex); }, {"ThetaCM", "Ex"});
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
    eff.Add("g0", "./Inputs/effs/g0_7Li_dp_sil.root", "effCM");
    //eff.Add("g1", "./Inputs/effs/g1_7Li_dp_sil.root", "effCM");
    // Draw to check is fine
    eff.Draw();

    // Set experiment info
    PhysUtils::Experiment exp {"../norm/7Li_norm.dat"};
    // And compute differential xs!
    Angular::DifferentialXS xs {&ivs, &fitter, &eff, &exp};
    xs.DoFor(peaks);
    xs.Write("./Outputs/");

    // Plot
    Angular::Comparator comp {"g.s", xs.Get("g0")};
    comp.Add("ADWA", "./Inputs/gs/21.gs");
    comp.Add("DA1p-Delaroche", "./Inputs/gs_DA1p_Delaroche/21.g0");
    comp.Add("Daehnik-Delaroche", "./Inputs/gs_Daehnik_Delaroche/21.g0");
    // comp.Add("DelarocheMe", "./Inputs/gs_Daehnik_Delaroche_myself/21.g0");
    comp.Fit();
    comp.Draw("", true);

    auto* c0 {new TCanvas {"c0", "(d,p) canvas"}};
    c0->DivideSquare(2);
    c0->cd(1);
    hEx->DrawClone();
    c0->cd(2);
    hCM->DrawClone("colz");
}
