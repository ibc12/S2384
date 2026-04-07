#include "ActCluster.h"
#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActLine.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActSRIM.h"
#include "ActTPCData.h"
#include "ActTPCParameters.h"
#include "ActVoxel.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"
#include <random>

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMath.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TSpline.h>
#include <TStyle.h>

#include <Math/Point3D.h>
#include <Math/Vector3D.h>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

#include "../../PrettyStyle.C"

////////////////////////////////////////////////////////
// Spline manually computing the Eloss from SRIM tables
////////////////////////////////////////////////////////

TSpline3* BuildSRIMsplineManual(ActPhysics::SRIM* srim, double range, const std::string& particleKey, double step = 0.5,
                                double sOffset = 0.0)
{
    TSpline3* sp = nullptr;
    std::vector<double> s_step_pts;
    std::vector<double> y_step_pts;

    double E = srim->EvalInverse(particleKey, range);

    bool ran_out = false;
    for(double r = 0; r < range; r += step)
    {
        double Epost = srim->Slow(particleKey, E, step);
        if(Epost < 0.0)
        {
            Epost = 0.0;
            ran_out = true;
        }
        double dE = E - Epost;
        E = Epost;

        if(dE <= 0)
        {
            if(ran_out)
                break;
            else
                continue;
        }

        double s = r + 0.5 * step + sOffset;
        s_step_pts.push_back(s);
        y_step_pts.push_back(dE);

        if(ran_out)
            break;
    }

    if(y_step_pts.back() > 0.0)
    {
        s_step_pts.push_back(range);
        y_step_pts.push_back(0.0);
    }

    int nSteps = (int)s_step_pts.size();
    if(nSteps < 3)
    {
        std::cout << "Not enough points to build spline (nSteps=" << nSteps << "). Need at least 3 points.\n";
        return nullptr;
    }
    sp = new TSpline3(("spSRIM_" + particleKey).c_str(), s_step_pts.data(), y_step_pts.data(), nSteps, "b1,e1", 0, 0);
    sp->SetNpx(3000);

    std::cout << "Spline max: " << sp->GetXmax() << "\n";
    std::cout << "Value near max: " << sp->Eval(range - 1) << "\n";

    return sp;
}

/////////////////////////////////////
// Spline computing the Eloss from SRIM tables, but using directly dE/dx to remove step size dependence
/////////////////////////////////////
TSpline3* buildSRIMspline(ActPhysics::SRIM* srim, double range, const std::string& particleKey, double step = 0.5,
                          double sOffset = 0.0)
{
    std::vector<double> s_pts;
    std::vector<double> y_pts;

    for(double s = 0; s < range; s += step)
    {
        double R = range - s; // remaining range

        if(R <= 0)
            break;

        double E = srim->EvalInverse(particleKey, R);
        double dEdx = srim->EvalStoppingPower(particleKey, E);

        if(dEdx < 0)
            dEdx = 0;

        s_pts.push_back(s + 0.5 * step + sOffset);
        y_pts.push_back(dEdx);
    }
    //  ensure endpoint goes to zero
    if(!y_pts.empty() && y_pts.back() > 0)
    {
        s_pts.push_back(range);
        y_pts.push_back(0.0);
    }

    if(s_pts.size() < 3)
    {
        std::cout << "Not enough points to build spline (nSteps=" << s_pts.size() << "). Need at least 3 points.\n";
        return nullptr;
    }

    auto* sp = new TSpline3(("spSRIM_" + particleKey).c_str(), s_pts.data(), y_pts.data(), s_pts.size());

    sp->SetNpx(3000);
    std::cout << "Spline max: " << sp->GetXmax() << "\n";
    std::cout << "Value near max: " << sp->Eval(range - 1) << "\n";
    return sp;
}

void CompareSPLINE_dEdxSRIM_dEdxManual()
{
    PrettyStyle(true, true);

    // Get dEdx from SRIM and from manual calculation (charge per voxel / voxel length) for all voxels in the light
    // cluster of L1 events; plot them against each other
        std::string beam {"7Li"};
        std::string target {"d"};
        std::string light {"p"};

        auto srim = new ActPhysics::SRIM();
        srim->ReadTable("light", "../../Calibrations/SRIM/1H_900mb_CF4_95-5.txt");

        auto spline_dEdxSRIM = buildSRIMspline(srim, 250, "light", 0.5);
        auto spline_dEdxManual = BuildSRIMsplineManual(srim, 250, "light", 0.5);

        auto* c = new TCanvas("c", "c", 800, 600);
        c->Divide(1, 2);
        c->cd(1);
        spline_dEdxSRIM->SetLineColor(kRed);
        spline_dEdxSRIM->Draw("AL");
        spline_dEdxManual->SetLineColor(kBlue);
        spline_dEdxManual->Draw("L SAME");
}