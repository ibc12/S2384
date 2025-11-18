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
#include "ActKinematics.h"

#include <string>
#include <vector>
#include "TStyle.h"
#include "TColor.h"

#include "../Histos.h"

void Ang(bool isLab = false)
{

    if(isLab)
        Angular::ToggleIsLab();

    ROOT::EnableImplicitMT();

    ROOT::RDataFrame df {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_7Li_d_d_filtered.root"};
    auto def {df.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == true; }, {"MergerData"})}; // only silicons

    // Book histograms
    auto hEx {def.Histo1D(S2384Fit::Exdd_7Li, "Ex")};
    ROOT::RDF::RResultPtr<TH2D> hKin {};
    if(isLab)
        hKin = def.Histo2D({"hKin", "Lab;#theta_{Lab};E_{Lab} [MeV]", 300, 0, 120, 300, 0, 60}, "fThetaLight", "EVertex");
    else
        hKin = def.Histo2D({"hCM", "CM;#theta_{CM};E [MeV]", 300, 0, 120, 300, 0, 60}, "ThetaCM", "EVertex");

    // Init intervals
    double thetaMin = isLab ? 55.0   : 34.0;
    double thetaMax = isLab ? 70.0 : 62.5;
    double thetaStep = isLab ? 2.5  : 2.50;
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
    for(const auto& peak : peaks)
    {
        TString inputPath = isLab ? TString::Format("Inputs/effs/%s_7Li_dd_sil_lab.root", peak.c_str())
                                        : TString::Format("Inputs/effs/%s_7Li_dd_sil.root", peak.c_str());
        eff.Add(peak, inputPath.Data(), isLab ? "effLab" : "effCM");
    }
    // Draw to check is fine
    eff.Draw();

    // Set experiment info
    PhysUtils::Experiment exp {4.126e19*25.6, 648000, 300}; // 2.08039e8
    // And compute differential xs!
    Angular::DifferentialXS xs {&ivs, &fitter, &eff, &exp};
    xs.DoFor(peaks);
    if(!isLab)
        xs.Write("./Outputs/");

    

    // Plot
    Angular::Comparator comp {"g.s", xs.Get("g0")};
    comp.Add("Haixia", "./Inputs/gsH/fort.201");
    comp.Add("Daehnick", "./Inputs/gsD/fort.201");
    comp.Add("DA1p", "./Inputs/gsDA1p/fort.201");
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

    auto* c0 {new TCanvas {"c0", "(d,d) canvas"}};
    c0->DivideSquare(2);
    c0->cd(1);
    hEx->DrawClone();
    c0->cd(2);
    hKin->DrawClone("colz");
}
