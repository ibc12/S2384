#include "ActMergerData.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TColor.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"

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

#include "../../PrettyStyle.C"
#include "../Histos.h"

void Ang()
{
    PrettyStyle(true, false);
    ROOT::EnableImplicitMT();

    ROOT::RDataFrame df {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_F_11Li_d_d_filtered.root"};
    // auto def {df.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == true; },
    //                     {"MergerData"})}; // only silicons, == false is for L1 events
    auto def {df.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; },
                        {"MergerData"})}; // false is for L1 events
    // Filter on some runs
    // auto def {
    //     df.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"})
    //         .Filter([](ActRoot::MergerData& m) { return m.fRun > 58; }, {"MergerData"})}; // false is for L1 events

    // Book histograms
    auto hEx {def.Histo1D(S2384Fit::Exdd, "Ex")};
    auto hCM {def.Histo2D({"hCM", "CM;#theta_{CM};E [MeV]", 300, 0, 120, 300, 0, 60}, "ThetaCM", "EVertex")};

    // Init intervals
    double thetaMin {12};
    double thetaMax {26};
    // double thetaMax {60};
    double thetaStep {3};
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
    // eff.Add("g0", "../../Simulation/Outputs/11Li/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0.root", "effCMside");
    eff.Add("g0",
            "../../Simulation/Outputs/11Li/test_triggerChange/"
            "2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_BothTrigger.root",
            "effCM");
    // eff.Add("g0",
    //         "../../Simulation/Outputs/11Li/test_triggerChange/"
    //         "2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_BeforeTrigger.root",
    //         "effCM");
    // eff.Add("g0",
    //         "../../Simulation/Outputs/11Li/test_triggerChange/"
    //         "2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_PostTrigger.root",
    //         "effCM");
    eff.Scale(0.931);
    // Draw to check is fine
    eff.Draw();

    // Set experiment info
    PhysUtils::Experiment exp {"../norm/11Li_norm_L1.dat"};
    // PhysUtils::Experiment exp {"../norm/norm_per_run_interval/11Li_norm_L1_runs59-122.dat"};
    // PhysUtils::Experiment exp {"../norm/norm_per_run_interval/11Li_norm_L1_runs19-58.dat"};
    // And compute differential xs!
    Angular::DifferentialXS xs {&ivs, &fitter, &eff, &exp};
    xs.DoFor(peaks);
    xs.Write("./Outputs/");

    // Plot
    Angular::Comparator comp {"g.s", xs.Get("g0")};
    // comp.Add("SFRESCO", "./Inputs/gsSFRESCO/fort.201");
    comp.Add("Haixia", "./Inputs/gsH/fort.201");
    // comp.Fit();
    comp.Add("Daehnick", "./Inputs/gsD/fort.201");
    // comp.Fit();
    comp.Add("DA1p", "./Inputs/gsDA1p/fort.201");
    comp.Add("full OP1", "../../../Tese/L11-Li12/full_OP1.dat");
    comp.Add("full OP2", "../../../Tese/L11-Li12/full_OP2.dat");
    comp.Add("no cont OP1", "../../../Tese/L11-Li12/no_cont_OP1.dat");
    comp.Add("no cont OP2", "../../../Tese/L11-Li12/no_cont_OP2.dat");
    comp.Fit();
    comp.Draw("", true);
    comp.DrawTheo();

    auto* c0 {new TCanvas {"c0", "(d,d) canvas"}};
    c0->DivideSquare(2);
    c0->cd(1);
    hEx->DrawClone();
    c0->cd(2);
    hCM->DrawClone("colz");
}
