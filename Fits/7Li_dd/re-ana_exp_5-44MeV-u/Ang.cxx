#include "ActKinematics.h"
#include "ActMergerData.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TColor.h"
#include "TGraphErrors.h"
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

#include "../../Histos.h"

// void Angular::Comparator::Flip(const std::string& name)
// {
//     auto graphs = this->GetTheoGraphs();
//     // verify graph exists
//     if(graphs.find(name) == graphs.end())
//         throw std::runtime_error("Comparator::Flip - Graph '" + name + "' not found");
//
//     // Get original graph and clone to avoid double-delete
//     auto gOrig = graphs.at(name);
//     auto gFlip = (TGraphErrors*) gOrig->Clone((name + "_flipped").c_str());
//
//     // Aply transformation θ' = 180 - θ
//     for(int i = 0; i < gFlip->GetN(); ++i)
//     {
//         double x, y;
//         gFlip->GetPoint(i, x, y);
//         gFlip->SetPoint(i, 180.0 - x, y);
//     }
//     // Sort from smallest to largest to be able to plot and replace it
//     gFlip->Sort();
//     this->Replace(name, gFlip);
// }

void Ang(bool isLab = false)
{

    if(isLab)
        Angular::ToggleIsLab();

    // paper data
    TGraphErrors* g_7Li_dd_gs_5_44MeV_u {new TGraphErrors("./7Li_dd_gs_5-44MeV-u.dat", "%lg %lg")};
    TGraphErrors* g_7Li_dp_gs_5_44MeV_u {new TGraphErrors("./7Li_dp_gs_5-44MeV-u.dat", "%lg %lg")};
    TGraphErrors* g_7Li_dp_1st_5_44MeV_u {new TGraphErrors("./7Li_dp_1st_5-44MeV-u.dat", "%lg %lg")};
    // My data
    TGraphErrors* g_7Li_dd_gs_7_5MeV_u {new TGraphErrors("../Outputs/xs/g0_xs.dat", "%lg %lg")};
    TGraphErrors* g_7Li_dp_gs_7_5MeV_u {new TGraphErrors("../../7Li_dp/Outputs/xs/g0_xs.dat", "%lg %lg")};
    TGraphErrors* g_7Li_dp_1st_7_5MeV_u {new TGraphErrors("../../7Li_dp/Outputs/xs/g1_xs.dat", "%lg %lg")};

    // Plot
    Angular::Comparator comp {"g.s paper", g_7Li_dd_gs_5_44MeV_u};
    comp.Add("DA1pcorr_MyEnergy", "../Inputs/gsDA1p_corr/fort.201");
    comp.Add("DA1pcorr_EnergyThisExp", "./Inputs/dd_DA1pcorr_5-44MeV-u/fort.201");
    comp.Fit();
    comp.Draw("", true);

    Angular::Comparator comp1 {"g.s mine", g_7Li_dd_gs_7_5MeV_u};
    comp1.Add("DA1pcorr_MyEnergy", "../Inputs/gsDA1p_corr/fort.201");
    comp1.Add("CH89_7_25_MeV_Ed", "./Inputs/dd_CH89_7_25MeV_Ed/fort.201");
    comp1.Add("CH89_14_5_MeV_Ed", "./Inputs/dd_CH89_14_5MeV_Ed/fort.201");
    comp1.Fit();
    comp1.Draw("", true);

    // Plot paper data and mine together
    // First for dd gs
    auto* c0 {new TCanvas {"c0", "(d,d') canvas"}};
    g_7Li_dd_gs_5_44MeV_u->SetMarkerStyle(20);
    g_7Li_dd_gs_5_44MeV_u->SetMarkerColor(TColor::GetColor("#d62728"));
    g_7Li_dd_gs_7_5MeV_u->SetMarkerStyle(21);
    g_7Li_dd_gs_7_5MeV_u->SetMarkerColor(TColor::GetColor("#9467bd"));
    g_7Li_dd_gs_5_44MeV_u->Draw("AP");
    g_7Li_dd_gs_7_5MeV_u->Draw("PSAME");
    // Now for dp gs
    auto* c1 {new TCanvas {"c1", "(d,p) gs canvas"}};
    g_7Li_dp_gs_5_44MeV_u->SetMarkerStyle(20);
    g_7Li_dp_gs_5_44MeV_u->SetMarkerColor(TColor::GetColor("#d62728"));
    g_7Li_dp_gs_7_5MeV_u->SetMarkerStyle(21);
    g_7Li_dp_gs_7_5MeV_u->SetMarkerColor(TColor::GetColor("#9467bd"));
    g_7Li_dp_gs_5_44MeV_u->Draw("AP");
    g_7Li_dp_gs_7_5MeV_u->Draw("PSAME");
    // Now for dp 1st
    auto* c2 {new TCanvas {"c2", "(d,p) 1st canvas"}};
    g_7Li_dp_1st_5_44MeV_u->SetMarkerStyle(20);
    g_7Li_dp_1st_5_44MeV_u->SetMarkerColor(TColor::GetColor("#d62728"));
    g_7Li_dp_1st_7_5MeV_u->SetMarkerStyle(21);
    g_7Li_dp_1st_7_5MeV_u->SetMarkerColor(TColor::GetColor("#9467bd"));
    g_7Li_dp_1st_5_44MeV_u->Draw("AP");
    g_7Li_dp_1st_7_5MeV_u->Draw("PSAME");
}
