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

#include <fstream>
#include <string>
#include <vector>

#include "../../PrettyStyle.C"
#include "../Histos.h"

void MergeSilL1()
{
    PrettyStyle();
    // Sil xs
    TGraphErrors* g_sil_xs {new TGraphErrors()};
    {
        std::ifstream in("./Outputs/xs/g0_xs.dat");

        double x, y, sy;
        while(in >> x >> y >> sy)
        {
            int n = g_sil_xs->GetN();
            g_sil_xs->SetPoint(n, x, y);
            g_sil_xs->SetPointError(n, 0., sy); // ex = 0, ey = sy
        }
    }

    // L1 xs
    TGraphErrors* g_L1_xs {new TGraphErrors()};
    {
        std::ifstream in("../dd_testL1/Outputs/xs/g0_xs.dat");

        double x, y, sy;
        while(in >> x >> y >> sy)
        {
            int n = g_L1_xs->GetN();
            g_L1_xs->SetPoint(n, x, y);
            g_L1_xs->SetPointError(n, 0., sy); // ex = 0, ey = sy
        }
    }

    // Create a TGraphErrors that merges both
    TGraphErrors* g_merged {new TGraphErrors()};
    int nSil = g_sil_xs->GetN();
    int nL1 = g_L1_xs->GetN();

    // Copy L1 points firstto have them sorted
    for(int i = 0; i < nL1; ++i)
    {
        double x, y, ex, ey;
        g_L1_xs->GetPoint(i, x, y);
        ex = g_L1_xs->GetErrorX(i);
        ey = g_L1_xs->GetErrorY(i);
        g_merged->SetPoint(g_merged->GetN(), x, y);
        g_merged->SetPointError(g_merged->GetN() - 1, ex, ey);
    }
    // Copy sil points
    for(int i = 0; i < nSil; ++i)
    {
        double x, y, ex, ey;
        g_sil_xs->GetPoint(i, x, y);
        ex = g_sil_xs->GetErrorX(i);
        ey = g_sil_xs->GetErrorY(i);
        g_merged->SetPoint(g_merged->GetN(), x, y);
        g_merged->SetPointError(g_merged->GetN() - 1, ex, ey);
    }

    Angular::Comparator comp {"Merged g.s", g_merged};
    // comp.Add("DA1pcorr", "./Inputs/gsDA1p_corr/fort.201");
    // comp.Fit();
    comp.Draw("", true);
    Angular::Comparator compSil {"Sil g.s", g_sil_xs};
    // compSil.Add("DA1pcorr", "./Inputs/gsDA1p_corr/fort.201");
    // compSil.Fit();
    compSil.Draw("", true);
    Angular::Comparator compL1 {"L1 g.s", g_L1_xs};
    // compL1.Add("DA1pcorr", "../7Li_dd/Inputs/gsDA1p_corr/fort.201");
    // compL1.Fit();
    compL1.Draw("", true);

    // Do the same for the 1st Ex
    // TGraphErrors* g_sil_xs_1stEx {new TGraphErrors()};
    // {
    //     std::ifstream in("./Outputs/xs/g1_xs.dat");
    //     double x, y, sy;
    //     while(in >> x >> y >> sy)
    //     {
    //         int n = g_sil_xs_1stEx->GetN();
    //         g_sil_xs_1stEx->SetPoint(n, x, y);
    //         g_sil_xs_1stEx->SetPointError(n, 0., sy); // ex = 0, ey = sy
    //     }
    // }
    // TGraphErrors* g_L1_xs_1stEx {new TGraphErrors()};
    // {
    //     std::ifstream in("../7Li_dd_testL1/Outputs/xs/g1_xs.dat");
    //     double x, y, sy;
    //     while(in >> x >> y >> sy)
    //     {
    //         int n = g_L1_xs_1stEx->GetN();
    //         g_L1_xs_1stEx->SetPoint(n, x, y);
    //         g_L1_xs_1stEx->SetPointError(n, 0., sy); // ex = 0, ey = sy
    //     }
    // }
    //
    // TGraphErrors* g_merged_1stEx {new TGraphErrors()};
    // int nSil_1stEx = g_sil_xs_1stEx->GetN();
    // int nL1_1stEx = g_L1_xs_1stEx->GetN();
    // // Copy L1 points firstto have them sorted
    // for(int i = 0; i < nL1_1stEx; ++i)
    // {
    //     double x, y, ex, ey;
    //     g_L1_xs_1stEx->GetPoint(i, x, y);
    //     ex = g_L1_xs_1stEx->GetErrorX(i);
    //     ey = g_L1_xs_1stEx->GetErrorY(i);
    //     g_merged_1stEx->SetPoint(g_merged_1stEx->GetN(), x, y);
    //     g_merged_1stEx->SetPointError(g_merged_1stEx->GetN() - 1, ex, ey);
    // }
    // // Copy sil points
    // for(int i = 0; i < nSil_1stEx; ++i)
    // {
    //     double x, y, ex, ey;
    //     g_sil_xs_1stEx->GetPoint(i, x, y);
    //     ex = g_sil_xs_1stEx->GetErrorX(i);
    //     ey = g_sil_xs_1stEx->GetErrorY(i);
    //     g_merged_1stEx->SetPoint(g_merged_1stEx->GetN(), x, y);
    //     g_merged_1stEx->SetPointError(g_merged_1stEx->GetN() - 1, ex, ey);
    // }
    //
    // Angular::Comparator comp_1stEx {"Merged 1st Ex", g_merged_1stEx};
    // comp_1stEx.Add("DA1pcorr BE2 deformation", "./Inputs/g1_DA1p_corr/fort.202");
    // comp_1stEx.Fit();
    // comp_1stEx.Draw("", true);
    // Angular::Comparator compSil_1stEx {"Sil 1st Ex", g_sil_xs_1stEx};
    // compSil_1stEx.Add("DA1pcorr BE2 deformation", "./Inputs/g1_DA1p_corr/fort.202");
    // compSil_1stEx.Fit();
    // compSil_1stEx.Draw("", true);
    // Angular::Comparator compL1_1stEx {"L1 1st Ex", g_L1_xs_1stEx};
    // compL1_1stEx.Add("DA1pcorr BE2 deformation", "../7Li_dd/Inputs/g1_DA1p_corr/fort.202");
    // compL1_1stEx.Fit();
    // compL1_1stEx.Draw("", true);
}