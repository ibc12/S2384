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

#include "../Histos.h"
#include "../../PrettyStyle.C"

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
        std::ifstream in("../7Li_dp_testL1/Outputs/xs/g0_xs.dat");

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

    // Do the same for the 1st excited state
    TGraphErrors* g_sil_xs_1st {new TGraphErrors()};
    {
        std::ifstream in("./Outputs/xs/g1_xs.dat");

        double x, y, sy;
        while(in >> x >> y >> sy)
        {
            int n = g_sil_xs_1st->GetN();
            g_sil_xs_1st->SetPoint(n, x, y);
            g_sil_xs_1st->SetPointError(n, 0., sy); // ex = 0, ey = sy
        }
    }

    TGraphErrors* g_L1_xs_1st {new TGraphErrors()};
    {
        std::ifstream in("../7Li_dp_testL1/Outputs/xs/g1_xs.dat");

        double x, y, sy;
        while(in >> x >> y >> sy)
        {
            int n = g_L1_xs_1st->GetN();
            g_L1_xs_1st->SetPoint(n, x, y);
            g_L1_xs_1st->SetPointError(n, 0., sy); // ex = 0, ey = sy
        }
    }
    TGraphErrors* g_merged_1st {new TGraphErrors()};
    int nSil_1st = g_sil_xs_1st->GetN();
    int nL1_1st = g_L1_xs_1st->GetN();
    // Copy L1 points firstto have them sorted
    for(int i = 0; i < nL1_1st; ++i)
    {
        double x, y, ex, ey;
        g_L1_xs_1st->GetPoint(i, x, y);
        ex = g_L1_xs_1st->GetErrorX(i);
        ey = g_L1_xs_1st->GetErrorY(i);
        g_merged_1st->SetPoint(g_merged_1st->GetN(), x, y);
        g_merged_1st->SetPointError(g_merged_1st->GetN() - 1, ex, ey);
    }
    // Copy sil points
    for(int i = 0; i < nSil_1st; ++i)
    {
        double x, y, ex, ey;
        g_sil_xs_1st->GetPoint(i, x, y);
        ex = g_sil_xs_1st->GetErrorX(i);
        ey = g_sil_xs_1st->GetErrorY(i);
        g_merged_1st->SetPoint(g_merged_1st->GetN(), x, y);
        g_merged_1st->SetPointError(g_merged_1st->GetN() - 1, ex, ey);
    }

    Angular::Comparator comp {"Merged g.s", g_merged};
    // comp.Add("DA1pcorr", "./Inputs/gs_DA1pcorr_Delaroche/21.g0");
    comp.Add("ADWA", "../7Li_dp/Inputs/gs_ADWA/fort.202");
    comp.Fit();
    comp.Draw("", true);
    Angular::Comparator compSil {"Sil g.s", g_sil_xs};
    compSil.Add("DA1pcorr", "./Inputs/gs_DA1pcorr_Delaroche/21.g0");
    compSil.Add("ADWA", "../7Li_dp/Inputs/gs_ADWA/fort.202");
    compSil.Fit();
    compSil.Draw("", true);
    Angular::Comparator compL1 {"L1 g.s", g_L1_xs};
    compL1.Add("DA1pcorr", "./Inputs/gs_DA1pcorr_Delaroche/21.g0");
    compL1.Add("ADWA", "../7Li_dp/Inputs/gs_ADWA/fort.202");
    compL1.Fit();
    compL1.Draw("", true);

    Angular::Comparator comp1 {"Merged 1st Ex", g_merged_1st};
    comp1.Add("DA1pcorr", "./Inputs/g1_DA1pcorr_Delaroche/21.g1");
    comp1.Add("ADWA", "../7Li_dp/Inputs/gs_ADWA/fort.203");
    comp1.Fit();
    comp1.Draw("", true);
    Angular::Comparator compSil1 {"Sil 1st Ex", g_sil_xs_1st};
    compSil1.Add("DA1pcorr", "./Inputs/g1_DA1pcorr_Delaroche/21.g1");
    compSil1.Add("ADWA", "../7Li_dp/Inputs/gs_ADWA/fort.203");
    compSil1.Fit();
    compSil1.Draw("", true);
    Angular::Comparator compL11 {"L1 1st Ex", g_L1_xs_1st};
    compL11.Add("DA1pcorr", "./Inputs/g1_DA1pcorr_Delaroche/21.g1");
    compL11.Add("ADWA", "../7Li_dp/Inputs/gs_ADWA/fort.203");
    compL11.Fit();
    compL11.Draw("", true);
}