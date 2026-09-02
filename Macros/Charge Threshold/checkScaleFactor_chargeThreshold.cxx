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

#include "../../Fits/Histos.h"
#include "../../PrettyStyle.C"

void checkScaleFactor_chargeThreshold(bool isLab = false)
{
    PrettyStyle(false);
    if(isLab)
        Angular::ToggleIsLab();

    ROOT::EnableImplicitMT();

    // ROOT::RDataFrame df {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_7Li_d_d_filtered.root"};
    ROOT::RDataFrame df {"Final_Tree", "./Inputs/tree_ex_F_7Li_elastic.root"}; // all elastic data, got from kinematic
                                                                               // plot Qtot vs ThetaLab
    auto def {df.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; },
                        {"MergerData"})}; // only silicons

    // Book histograms
    auto hEx {def.Histo1D(S2384Fit::Exdd_7Li, "Ex")};
    ROOT::RDF::RResultPtr<TH2D> hKin {};
    if(isLab)
        hKin =
            def.Histo2D({"hKin", "Lab;#theta_{Lab};E_{Lab} [MeV]", 300, 0, 120, 300, 0, 60}, "fThetaLight", "EVertex");
    else
        hKin = def.Histo2D({"hCM", "CM;#theta_{CM};E [MeV]", 300, 0, 120, 300, 0, 60}, "ThetaCM", "EVertex");

    // Init intervals
    double thetaMin = isLab ? 70.0 : 16;
    double thetaMax = isLab ? 80.0 : 26;
    double thetaStep = isLab ? 2.5 : 2.50;
    Angular::Intervals ivs {thetaMin, thetaMax, S2384Fit::Exdd_7Li, thetaStep, 0};
    if(isLab)
        def.Foreach([&](float thetalab, double ex) { ivs.Fill(thetalab, ex); }, {"fThetaLight", "Ex"});
    else
        def.Foreach([&](double thetacm, double ex) { ivs.Fill(thetacm, ex); }, {"ThetaCM", "Ex"});
    ivs.Draw();

    // Init fitter
    Angular::Fitter fitter {&ivs};
    // fitter.SetAllowFreeMean(true);
    // fitter.SetFreeMeanRange(0.1);
    fitter.Configure("../../Fits/7Li_dd_testL1/Outputs/fit.root");
    fitter.Run();
    fitter.Draw();
    fitter.DrawCounts();

    // Interface
    Fitters::Interface inter;
    inter.Read("../../Fits/7Li_dd_testL1/Outputs/interface.root");

    auto peaks {inter.GetKeys()};

    // Efficiency
    Interpolators::Efficiency eff;
    // for(const auto& peak : peaks)
    //{
    //     TString inputPath = isLab ? TString::Format("Inputs/effs/%s_7Li_dd_sil_lab.root", peak.c_str())
    //                                     : TString::Format("Inputs/effs/%s_7Li_dd_sil.root", peak.c_str());
    //     eff.Add(peak, inputPath.Data(), isLab ? "effLab" : "effCM");
    // }
    // eff.Add("g0", "../../Simulation/Outputs/7Li/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1.root", isLab ? "effLab" :
    // "effCM"); eff.Add("g1", "../../Simulation/Outputs/7Li/2H_2H_TRIUMF_Eex_0.477_nPS_0_pPS_0_L1.root", isLab ?
    // "effLab" : "effCM");
    // eff.Add("g0",
    // "../../Simulation/Outputs/7Li/test_charge_threshold/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_1e6Thresh.root", isLab
    // ? "effLab" : "effCM"); eff.Add("g1",
    // "../../Simulation/Outputs/7Li/test_charge_threshold/2H_2H_TRIUMF_Eex_0.477_nPS_0_pPS_0_L1_1e6Thresh.root", isLab
    // ? "effLab" : "effCM");
    eff.Add("g0",
            "../../Simulation/Outputs/7Li/test_charge_threshold/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_1e5Thresh.root",
            isLab ? "effLab" : "effCM");
    eff.Scale(0.9);
    eff.Add("g1",
            "../../Simulation/Outputs/7Li/test_charge_threshold/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_1e5Thresh.root",
            isLab ? "effLab" : "effCM");
    // Draw to check is fine
    eff.Draw();

    // Set experiment info
    PhysUtils::Experiment exp {"../../Fits/norm/7Li_norm.dat"};
    // And compute differential xs!
    Angular::DifferentialXS xs {&ivs, &fitter, &eff, &exp};
    xs.DoFor(peaks);

    // Plot
    Angular::Comparator comp {"g.s", xs.Get("g0")};
    comp.Add("Haixia", "../../Fits/7Li_dd/Inputs/gsH/fort.201");
    comp.Add("Daehnick", "../../Fits/7Li_dd/Inputs/gsD/fort.201");
    comp.Add("DA1p", "../../Fits/7Li_dd/Inputs/gsDA1p/fort.201");
    comp.Add("DA1pcorr", "../../Fits/7Li_dd/Inputs/gsDA1p_corr/fort.201");
    // comp.Add("ADWA", "../7Li_dp/Inputs/gs_ADWA/fort.201");
    Angular::Comparator comp1 {"1st Ex", xs.Get("g1")};
    comp1.Add("DA1p BE2 deformation", "../../Fits/7Li_dd/Inputs/g1_DA1p/fort.202");
    comp1.Add("DA1pcorr BE2 deformation", "../../Fits/7Li_dd/Inputs/g1_DA1p_corr/fort.202");
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
    std::cout << "Scale factor for DA1pcorr: " << comp.GetSF("DA1pcorr") << std::endl;

    // ---- Scan de charge threshold vs Scale Factor ----
    std::vector<std::string> threshLabels {"1e4", "2e4", "3e4", "4e4", "5e4", "6e4", "7e4", "8e4", "9e4", "1e5"};
    std::vector<double> threshValues {1e4, 2e4, 3e4, 4e4, 5e4, 6e4, 7e4, 8e4, 9e4, 1e5};
    std::vector<double> sfValues;

    for(const auto& thresh : threshLabels)
    {
        // Reconstruir eficiencia para este threshold
        Interpolators::Efficiency effLoop;
        TString path0 = TString::Format(
            "../../Simulation/Outputs/7Li/test_charge_threshold/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_%sThresh.root",
            thresh.c_str());
        effLoop.Add("g0", path0.Data(), isLab ? "effLab" : "effCM");
        effLoop.Scale(0.9);
        effLoop.Add("g1", path0.Data(), isLab ? "effLab" : "effCM");

        // Recalcular xs con esta eficiencia
        Angular::DifferentialXS xsLoop {&ivs, &fitter, &effLoop, &exp};
        xsLoop.DoFor(peaks);

        // Rehacer el comparador solo para g.s. (g0) y sacar el SF de DA1pcorr
        Angular::Comparator compLoop {"g.s", xsLoop.Get("g0")};
        compLoop.Add("Haixia", "../../Fits/7Li_dd/Inputs/gsH/fort.201");
        compLoop.Add("Daehnick", "../../Fits/7Li_dd/Inputs/gsD/fort.201");
        compLoop.Add("DA1p", "../../Fits/7Li_dd/Inputs/gsDA1p/fort.201");
        compLoop.Add("DA1pcorr", "../../Fits/7Li_dd/Inputs/gsDA1p_corr/fort.201");
        compLoop.Fit();

        double sf = compLoop.GetSF("DA1pcorr");
        sfValues.push_back(sf);

        std::cout << "Threshold = " << thresh << "  ->  SF(DA1pcorr) = " << sf << std::endl;
    }

    // Gráfica SF vs charge threshold
    auto* gSF {new TGraph((int)threshValues.size(), threshValues.data(), sfValues.data())};
    gSF->SetTitle("Scale factor vs charge threshold;Charge threshold;Scale factor (DA1pcorr)");
    gSF->SetMarkerStyle(20);
    gSF->SetMarkerColor(kBlue + 1);
    gSF->SetLineColor(kBlue + 1);
    gSF->SetLineWidth(2);

    auto* cSF {new TCanvas {"cSF", "SF vs charge threshold"}};
    gSF->Draw("APL");

    // comp.ScaleToExp("DA1pcorr", &exp, xs.Get("g0"), eff.GetTEfficiency("g0"), 0.080);
    comp.QuotientPerPoint();

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