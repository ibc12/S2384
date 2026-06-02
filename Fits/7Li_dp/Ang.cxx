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

void Ang(bool isLab = false)
{
    PrettyStyle(false);
    if(isLab)
        Angular::ToggleIsLab();

    ROOT::EnableImplicitMT();

    ROOT::RDataFrame df {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_7Li_d_p_filtered.root"};
    auto def {df.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == true; },
                        {"MergerData"})}; // only silicons, == false is for L1 events
    // Book histograms
    auto hEx {def.Histo1D(S2384Fit::Exdp_7Li, "Ex")};
    auto hCM {def.Histo2D({"hCM", "CM;#theta_{CM};E [MeV]", 300, 0, 120, 300, 0, 60}, "ThetaCM", "EVertex")};

    // Init intervals
    double thetaMin {32};
    double thetaMax {63};
    double thetaStep {5};
    Angular::Intervals ivs {thetaMin, thetaMax, S2384Fit::Exdp_7Li, thetaStep, 0};
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
    eff.Add("g0", "../../Simulation/Outputs/7Li/2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0.root", isLab ? "effLab" : "effCM");
    eff.Add("g1", "../../Simulation/Outputs/7Li/2H_1H_TRIUMF_Eex_0.981_nPS_0_pPS_0.root", isLab ? "effLab" : "effCM");
    // eff.Add("g2", "../../Simulation/Outputs/7Li/2H_1H_TRIUMF_Eex_2.255_nPS_0_pPS_0.root", isLab ? "effLab" : "effCM");
    // eff.Add("v0", "../../Simulation/Outputs/7Li/2H_1H_TRIUMF_Eex_3.210_nPS_0_pPS_0.root", isLab ? "effLab" : "effCM");
    // eff.Add("v1", "../../Simulation/Outputs/7Li/2H_1H_TRIUMF_Eex_5.400_nPS_0_pPS_0.root", isLab ? "effLab" : "effCM");
    // eff.Add("v2", "../../Simulation/Outputs/7Li/2H_1H_TRIUMF_Eex_6.100_nPS_0_pPS_0.root", isLab ? "effLab" : "effCM");
    // eff.Add("v3", "../../Simulation/Outputs/7Li/2H_1H_TRIUMF_Eex_6.530_nPS_0_pPS_0.root", isLab ? "effLab" : "effCM");
    // eff.Add("v4", "../../Simulation/Outputs/7Li/2H_1H_TRIUMF_Eex_7.100_nPS_0_pPS_0.root", isLab ? "effLab" : "effCM");
    // eff.Add("ps0", "../../Simulation/Outputs/7Li/2H_1H_TRIUMF_Eex_0.000_nPS_1_pPS_0.root", isLab ? "effLab" : "effCM");
    // eff.Add("g1", "./Inputs/effs/g1_7Li_dp_sil.root", "effCM");
    //  Draw to check is fine
    eff.Draw();

    // Set experiment info
    PhysUtils::Experiment exp {"../norm/7Li_norm.dat"};
    // And compute differential xs!
    Angular::DifferentialXS xs {&ivs, &fitter, &eff, &exp};
    xs.DoFor(peaks);
    xs.Write("./Outputs/");

    // Plot
    Angular::Comparator comp {"g.s", xs.Get("g0")};
    comp.Add("ADWA", "./Inputs/gs_ADWA/fort.202");
    comp.Add("ADWA twofnr", "./Inputs/gs_ADWA/21.gs");
    comp.Add("ADWA finite range", "./Inputs/gs_ADWA_FiniteRange/fort.202");
    comp.Add("DA1p-Delaroche", "./Inputs/gs_DA1p_Delaroche/21.g0");
    comp.Add("Daehnik-Delaroche", "./Inputs/gs_Daehnik_Delaroche/21.g0");
    comp.Add("DA1pcorr-Delaroche", "./Inputs/gs_DA1pcorr_Delaroche/21.g0");
    comp.Add("DA1pcorr-Delaroche tweakPar", "./Inputs/gs_DA1pcorr_Delaroche_tweakPar/21.g0");
    comp.Add("DA1pcorr-JLMmicroscopic (fermi)", "./Inputs/gs_DA1pcorr_JLMmicroscopic/21.g0");
    comp.Add("DA1pcorr-JLMmicroscopic (oscillator)", "./Inputs/gs_DA1pcorr_JLMmicroscopic_other/21.g0");
    comp.Add("DA1pcorr-Perey", "./Inputs/gs_DA1pcorr_Perey/21.g0");
    comp.Add("DA1pcorr-CH89", "./Inputs/gs_DA1pcorr_CH89/21.g0");
    comp.Add("DA1pcorr-Becheti", "./Inputs/gs_DA1pcorr_Bechetti/21.g0");
    comp.Add("DA1pcorr-Menet", "./Inputs/gs_DA1pcorr_Menet/21.g0");
    comp.Add("OMP paper 5.44 MeV/u", "./Inputs/gs_DWBA_paper/21.g0");
    comp.Add("ADWA-CH89", "./Inputs/gs_ADWA_CH89/21.g0");
    comp.Add("ADWA-Perey", "./Inputs/gs_ADWA_Perey/21.g0");
    comp.Fit();
    comp.Draw("", true);
    comp.DrawTheo();
    comp.DrawSFfromIntegral();

    auto* c0 {new TCanvas {"c0", "(d,p) canvas"}};
    c0->DivideSquare(2);
    c0->cd(1);
    hEx->DrawClone();
    c0->cd(2);
    hCM->DrawClone("colz");

    double chi2Intervals {};
    for(const auto& result : fitter.GetTFitResults())
    {
        // Process fit results
        chi2Intervals += result.Chi2();
    }
    std::cout << "Total chi2: " << chi2Intervals / fitter.GetTFitResult(0).Ndf() << "\n";
}
