#include "ActKinematics.h"
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

void Ang_lDependent(bool isLab = false)
{
    // PrettyStyle(false);
    if(isLab)
        Angular::ToggleIsLab();

    ROOT::EnableImplicitMT();

    ROOT::RDataFrame df {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_F_7Li_d_p_filtered.root"};
    // ROOT::RDataFrame df {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_F_7Li_d_p_filtered.root"};
    auto def {df.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; },
                        {"MergerData"})}; // only silicons

    // Book histograms
    auto hEx {def.Histo1D(S2384Fit::Exdp_7Li, "Ex")};

    ROOT::RDF::RResultPtr<TH2D> hKin {};
    if(isLab)
        hKin =
            def.Histo2D({"hKin", "Lab;#theta_{Lab};E_{Lab} [MeV]", 300, 0, 120, 300, 0, 60}, "fThetaLight", "EVertex");
    else
        hKin = def.Histo2D({"hCM", "CM;#theta_{CM};E [MeV]", 300, 0, 120, 300, 0, 60}, "ThetaCM", "EVertex");

    // Init intervals
    double thetaMin = isLab ? 70.0 : 8;
    double thetaMax = isLab ? 80.0 : 40;
    double thetaStep = isLab ? 2.5 : 2.0;
    Angular::Intervals ivs {thetaMin, thetaMax, S2384Fit::Exdp_7Li, thetaStep, 2};
    if(isLab)
        def.Foreach([&](float thetalab, double ex) { ivs.Fill(thetalab, ex); }, {"fThetaLight", "Ex"});
    else
        def.Foreach([&](double thetacm, double ex) { ivs.Fill(thetacm, ex); }, {"ThetaCM", "Ex"});

    // Phase space
    ROOT::RDataFrame phase {"SimulationTTree", "../../Simulation/Outputs/7Li/test_rpCut/"
                                               "2H_1H_TRIUMF_Eex_0.000_nPS_1_pPS_0_L1_240rpCut_paralelization.root"};
    auto hPS {phase.Histo1D(S2384Fit::Exdp_7Li, "Eex", "weight")};
    hPS->SetNameTitle("hPS", "1n PS");
    phase.Foreach([&](double thetacm, double ex, double weight) { ivs.FillPS(0, thetacm, ex, weight); },
                  {"theta3CM", "Eex", "weight"});
    ivs.Draw();
    // Histogram of contamination
    ROOT::RDataFrame phase2 {"Final_Tree", "../../Macros/SimulationExp Comparison/Outputs/tree_ex_F_7Li_d_as_p.root"};
    auto hPS2 {phase2.Histo1D(S2384Fit::Exdp_7Li, "Ex_d_as_p")};
    hPS2->SetNameTitle("hPS2", "1n PS2");
    phase2.Foreach([&](double thetacm, double ex) { ivs.FillPS(1, thetacm, ex, 1); }, {"ThetaCM", "Ex_d_as_p"});
    ivs.TreatPS(0, 0.2);
    ivs.Draw();

    // Init fitter
    Angular::Fitter fitter {&ivs};
    // fitter.SetAllowFreeMean(true);
    // fitter.SetFreeMeanRange(0.1);
    fitter.Configure("./Outputs/fit_lDependent.root");
    double R {(std::pow(7, 1. / 3.) + std::pow(1, 1. / 3.)) * 1.25}; // interaction radius in fm, with r0 = 1.25 fm
    double mu {7 * 1. / (7 + 1) * 931.5};                            // reduced mass in MeV/c^2
    fitter.ApplyLambdaToModels(
        [mu, R](Fitters::Model& m)
        {
            m.AddBWL(0, 1, 2.03262, mu, R);
            m.AddBWL(1, 1, 2.03262, mu, R);
            m.AddBWL(2, 1, 2.03262, mu, R);
            m.AddBWL(3, 2, 2.03262, mu, R);
            m.AddBWL(4, 2, 2.03262, mu, R);
            m.AddBWL(5, 1, 2.03262, mu, R);
            m.AddBWL(6, 2, 2.03262, mu, R);
            m.AddBWL(7, 2, 2.03262, mu, R);
            m.AddBWL(8, 2, 2.03262, mu, R);
            m.AddBWL(9, 1, 2.03262, mu, R);
        });
    fitter.Run();
    fitter.Draw();
    fitter.DrawCounts();

    // Interface
    Fitters::Interface inter;
    inter.Read("./Outputs/interface_lDependent.root");

    auto peaks {inter.GetKeys()};

    // Efficiency
    Interpolators::Efficiency eff;
    // for(const auto& peak : peaks)
    //{
    //     TString inputPath = isLab ? TString::Format("Inputs/effs/%s_7Li_dd_sil_lab.root", peak.c_str())
    //                                     : TString::Format("Inputs/effs/%s_7Li_dd_sil.root", peak.c_str());
    //     eff.Add(peak, inputPath.Data(), isLab ? "effLab" : "effCM");
    // }
    // eff.Add("g0", "../../Simulation/Outputs/7Li/2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1.root",
    //         isLab ? "effLab" : "effCM");
    // eff.Add("g0",
    //         "../../Simulation/Outputs/7Li/test_ang_straggling_L1/"
    //         "2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_3-5AngStr_14pads.root",
    //         isLab ? "effLab" : "effCM");
    // eff.Add("g1",
    //         "../../Simulation/Outputs/7Li/test_ang_straggling_L1/"
    //         "2H_1H_TRIUMF_Eex_0.981_nPS_0_pPS_0_L1_3-5AngStr_14pads.root",
    //         isLab ? "effLab" : "effCM");
    eff.Add("g0",
            "../../Simulation/Outputs/7Li/test_rpCut/"
            "2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_240rpCut_paralelization.root",
            isLab ? "effLab" : "effCM");
    eff.Add("g1",
            "../../Simulation/Outputs/7Li/test_rpCut/"
            "2H_1H_TRIUMF_Eex_0.981_nPS_0_pPS_0_L1_240rpCut_paralelization.root",
            isLab ? "effLab" : "effCM");
    eff.Add("v0",
            "../../Simulation/Outputs/7Li/test_rpCut/"
            "2H_1H_TRIUMF_Eex_2.255_nPS_0_pPS_0_L1_240rpCut_paralelization.root",
            isLab ? "effLab" : "effCM");
    eff.Add("v1",
            "../../Simulation/Outputs/7Li/test_rpCut/"
            "2H_1H_TRIUMF_Eex_3.210_nPS_0_pPS_0_L1_240rpCut_paralelization.root",
            isLab ? "effLab" : "effCM");
    eff.Add("v2",
            "../../Simulation/Outputs/7Li/test_rpCut/"
            "2H_1H_TRIUMF_Eex_4.100_nPS_0_pPS_0_L1_240rpCut_paralelization.root",
            isLab ? "effLab" : "effCM");
    eff.Add("v3",
            "../../Simulation/Outputs/7Li/test_rpCut/"
            "2H_1H_TRIUMF_Eex_5.400_nPS_0_pPS_0_L1_240rpCut_paralelization.root",
            isLab ? "effLab" : "effCM");
    eff.Add("v4",
            "../../Simulation/Outputs/7Li/test_rpCut/"
            "2H_1H_TRIUMF_Eex_6.100_nPS_0_pPS_0_L1_240rpCut_paralelization.root",
            isLab ? "effLab" : "effCM");
    eff.Add("v5",
            "../../Simulation/Outputs/7Li/test_rpCut/"
            "2H_1H_TRIUMF_Eex_6.530_nPS_0_pPS_0_L1_240rpCut_paralelization.root",
            isLab ? "effLab" : "effCM");
    eff.Add("v6",
            "../../Simulation/Outputs/7Li/test_rpCut/"
            "2H_1H_TRIUMF_Eex_7.100_nPS_0_pPS_0_L1_240rpCut_paralelization.root",
            isLab ? "effLab" : "effCM");
    eff.Add("v7",
            "../../Simulation/Outputs/7Li/test_rpCut/"
            "2H_1H_TRIUMF_Eex_7.100_nPS_0_pPS_0_L1_240rpCut_paralelization.root",
            isLab ? "effLab" : "effCM");
    eff.Add("v8",
            "../../Simulation/Outputs/7Li/test_rpCut/"
            "2H_1H_TRIUMF_Eex_7.100_nPS_0_pPS_0_L1_240rpCut_paralelization.root",
            isLab ? "effLab" : "effCM");
    eff.Add("v9",
            "../../Simulation/Outputs/7Li/test_rpCut/"
            "2H_1H_TRIUMF_Eex_2.255_nPS_0_pPS_0_L1_240rpCut_paralelization.root",
            isLab ? "effLab" : "effCM");
    eff.Scale(0.931);
    // Draw to check is fine
    eff.Draw();

    // Set experiment info
    PhysUtils::Experiment exp {"../norm/7Li_norm_L1.dat"};
    // And compute differential xs!
    Angular::DifferentialXS xs {&ivs, &fitter, &eff, &exp};
    xs.DoFor(peaks);
    xs.TrimX("g0", 16, false);
    xs.TrimX("g0", 8, true);
    xs.TrimX("g1", 19, false);
    xs.TrimX("g1", 10, true);
    xs.TrimX("v0", 23, false);
    xs.TrimX("v0", 13, true);
    xs.TrimX("v1", 25, false);
    xs.TrimX("v1", 18, true);
    xs.TrimX("v2", 28, false);
    xs.TrimX("v2", 18, true);
    xs.TrimX("v3", 31, false);
    xs.TrimX("v3", 23, true);
    xs.TrimX("v4", 32, false);
    xs.TrimX("v4", 25, true);
    xs.TrimX("v5", 33, false);
    xs.TrimX("v5", 25, true);
    xs.TrimX("v6", 35, false);
    xs.TrimX("v6", 28, true);
    if(!isLab)
        xs.Write("./Outputs/", "xs_lDependent");

    // Plot
    Angular::Comparator comp {"g.s", xs.Get("g0")};
    comp.Add("DA1p-Delaroche", "./Inputs/gs_DA1pcorr_Delaroche/21.g0");
    comp.Add("ADWA", "../7Li_dp/Inputs/gs_ADWA/fort.202");
    comp.Add("ADWA DWUCK", "../FranckCalculations/DWUCK4_personal_calculations/good_calculations/gs/xs_dw4.n");
    Angular::Comparator comp1 {"1st Ex", xs.Get("g1")};
    comp1.Add("DA1p-Delaroche", "./Inputs/g1_DA1pcorr_Delaroche/21.g1");
    comp1.Add("ADWA", "../7Li_dp/Inputs/gs_ADWA/fort.203");
    comp1.Add("ADWA DWUCK", "../FranckCalculations/DWUCK4_personal_calculations/good_calculations/g1/xs_dw4.n");
    Angular::Comparator comp2 {"2nd Ex", xs.Get("v0")};
    comp2.Add("ADWA", "../7Li_dp/Inputs/gs_ADWA/fort.204");
    comp2.Add("ADWA DWUCK", "../FranckCalculations/DWUCK4_personal_calculations/good_calculations/v0/xs_dw4.n");
    Angular::Comparator comp3 {"3rd Ex", xs.Get("v1")};
    comp3.Add("ADWA DWUCK", "../FranckCalculations/DWUCK4_personal_calculations/good_calculations/v1/xs_dw4.n");
    // Plot
    if(isLab)
    {
        ActPhysics::Kinematics kin {"7Li(d,d)@51"};
        for(const auto& peak : peaks)
        {
            auto theo {comp.GetTheoGraphs()};
            for(const auto& [name, gtheo] : theo)
            {
                auto trans {kin.TransfromCMCrossSectionToLab(gtheo)};
                comp.Replace(name, (TGraphErrors*)trans->Clone());
                delete trans;
            }
        }
    }
    comp.Fit();
    comp.Draw("", true);
    comp.DrawTheo();
    comp1.Fit();
    comp1.Draw("", true);
    comp1.DrawTheo();
    comp2.Fit();
    comp2.Draw("", true);
    comp2.DrawTheo();
    comp3.Fit();
    comp3.Draw("", true);
    comp3.DrawTheo();
    auto* c0 {new TCanvas {"c0", "(d,d) canvas"}};
    c0->DivideSquare(2);
    c0->cd(1);
    hEx->DrawClone();
    c0->cd(2);
    hKin->DrawClone("colz");

    double chi2Intervals {};
    for(const auto& result : fitter.GetTFitResults())
    {
        // Process fit results
        chi2Intervals += result.Chi2();
    }
    std::cout << "Total chi2: " << chi2Intervals / fitter.GetTFitResult(0).Ndf() << "\n";
}
