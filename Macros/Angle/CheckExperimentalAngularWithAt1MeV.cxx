#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActTPCData.h"

#include <ROOT/RDataFrame.hxx>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TTree.h"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

void CheckExperimentalAngularWithAt1MeV()
{
    // Get dataframes for 7Li dp and dd
    ROOT::EnableImplicitMT();

    ROOT::RDataFrame df_dp {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_F_7Li_d_p_filtered.root"};
    auto def_dp {
        df_dp.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"})}; // only L1

    ROOT::RDataFrame df_dd {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_F_7Li_d_d_filtered.root"};
    auto def_dd {
        df_dd.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"})}; // only L1

    // Cut both datafrmames to Ex < 0.5 MeV to have only gs
    auto def_dp_gs {def_dp.Filter([](double ex) { return ex < 0.5; }, {"Ex"})};
    auto def_dd_gs {def_dd.Filter([](double ex) { return ex < 0.5; }, {"Ex"})};

    // Get a fraction of the energy region to see the angular dispersion
    double Emin = 1.4;
    double Emax = 1.6;
    auto def_dp_1MeV {def_dp_gs.Filter([=](double EVertex) { return EVertex > Emin && EVertex < Emax; }, {"EVertex"})};
    auto def_dd_1MeV {def_dd_gs.Filter([=](double EVertex) { return EVertex > Emin && EVertex < Emax; }, {"EVertex"})};

    // Plot m.fThetaLight for both dataframes
    auto hTheta_dp_gs {def_dp_gs.Histo1D({"hTheta_dp_gs", "7Li(d,p) gs;#theta_{Lab} [deg];Counts", 100, 70, 180},
                                         "MergerData.fThetaLight")};
    auto hTheta_dd_gs {def_dd_gs.Histo1D({"hTheta_dd_gs", "7Li(d,d) gs;#theta_{Lab} [deg];Counts", 100, 70, 180},
                                         "MergerData.fThetaLight")};
    auto hTheta_dp_1MeV {def_dp_1MeV.Histo1D(
        {"hTheta_dp_1MeV", "7Li(d,p) 1.4-1.6 MeV;#theta_{Lab} [deg];Counts", 100, 70, 180}, "MergerData.fThetaLight")};
    auto hTheta_dd_1MeV {def_dd_1MeV.Histo1D(
        {"hTheta_dd_1MeV", "7Li(d,d) 1.4-1.6 MeV;#theta_{Lab} [deg];Counts", 100, 70, 180}, "MergerData.fThetaLight")};

    // Now the kinematic EVertex vs theta for both dataframes
    auto hKin_dp_gs {
        def_dp_gs.Histo2D({"hKin_dp_gs", "7Li(d,p) gs;#theta_{Lab} [deg];E_{Lab} [MeV]", 100, 70, 180, 100, 0, 5},
                          "MergerData.fThetaLight", "EVertex")};
    auto hKin_dd_gs {
        def_dd_gs.Histo2D({"hKin_dd_gs", "7Li(d,d) gs;#theta_{Lab} [deg];E_{Lab} [MeV]", 100, 70, 180, 100, 0, 5},
                          "MergerData.fThetaLight", "EVertex")};
    auto hKin_dp_1MeV {def_dp_1MeV.Histo2D(
        {"hKin_dp_1MeV", "7Li(d,p) 1.4-1.6 MeV;#theta_{Lab} [deg];E_{Lab} [MeV]", 100, 70, 180, 100, 0, 5},
        "MergerData.fThetaLight", "EVertex")};
    auto hKin_dd_1MeV {def_dd_1MeV.Histo2D(
        {"hKin_dd_1MeV", "7Li(d,d) 1.4-1.6 MeV;#theta_{Lab} [deg];E_{Lab} [MeV]", 100, 70, 180, 100, 0, 5},
        "MergerData.fThetaLight", "EVertex")};

    // Now draw the histograms in a canvas
    TCanvas* c1 = new TCanvas("c1", "Angular distributions", 1200, 800);
    c1->Divide(2, 2);
    c1->cd(1);
    hTheta_dp_gs->DrawClone();
    c1->cd(2);
    hTheta_dd_gs->DrawClone();
    c1->cd(3);
    hTheta_dp_1MeV->DrawClone();
    c1->cd(4);
    hTheta_dd_1MeV->DrawClone();

    TCanvas* c2 = new TCanvas("c2", "Kinematic distributions", 1200, 800);
    c2->Divide(2, 2);
    c2->cd(1);
    hKin_dp_gs->DrawClone("colz");
    c2->cd(2);
    hKin_dd_gs->DrawClone("colz");
    c2->cd(3);
    hKin_dp_1MeV->DrawClone("colz");
    c2->cd(4);
    hKin_dd_1MeV->DrawClone("colz");
}
