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
#include <fstream>

#include "../Histos.h"

void MergeSilL1()
{
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
        std::ifstream in("../7Li_dd_testL1/Outputs/xs/g0_xs.dat");

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
    comp.Add("DA1pcorr", "./Inputs/gsDA1p_corr/fort.201");
    comp.Fit();
    comp.Draw("", true);
    Angular::Comparator compSil {"Sil g.s", g_sil_xs};
    compSil.Add("DA1pcorr", "./Inputs/gsDA1p_corr/fort.201");
    compSil.Fit();
    compSil.Draw("", true);
    Angular::Comparator compL1 {"L1 g.s", g_L1_xs};
    compL1.Add("DA1pcorr", "../7Li_dd/Inputs/gsDA1p_corr/fort.201");
    compL1.Fit();
    compL1.Draw("", true);
}