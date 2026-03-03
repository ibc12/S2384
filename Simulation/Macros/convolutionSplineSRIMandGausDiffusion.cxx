#include "ActCluster.h"
#include "ActLine.h"
#include "ActSRIM.h"
#include "ActTPCParameters.h"
#include "ActVoxel.h"

#include <random>

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TKey.h>
#include <TLegend.h>
#include <TMath.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TSpline.h>
#include <TStyle.h>

#include <Math/Point3D.h>
#include <Math/Vector3D.h>
#include <cmath>
#include <map>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

#include "../../PrettyStyle.C"

// ------------------------------------------------------------
// Building spline from SRIM sampling by steps
// ------------------------------------------------------------
TSpline3* BuildSRIMspline(ActPhysics::SRIM* srim, double range, const std::string& particleKey, double step = 0.5)
{
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

        std::cout << "r = " << r << ", E = " << E << ", dE = " << dE << "\n";

        if(dE <= 0)
        {
            if(ran_out)
            {
                s_step_pts.push_back(r);
                y_step_pts.push_back(0.0);
                break;
            }
            else
                continue;
        }

        double s = r + 0.5 * step;
        s_step_pts.push_back(s);
        y_step_pts.push_back(dE);

        if(ran_out)
        {
            s_step_pts.push_back(r + step);
            y_step_pts.push_back(0.0);
            break;
        }
    }

    std::cout << "Final values: s = " << s_step_pts.back() << ", dE = " << y_step_pts.back() << "\n";

    if(y_step_pts.back() > 0.0)
    {
        s_step_pts.push_back(range);
        y_step_pts.push_back(0.0);
    }

    // -------------------------------------------------
    // 3) Build spline from reduced set of significant points
    // -------------------------------------------------

    double slope_end = (0.0 - y_step_pts[y_step_pts.size() - 2]) / (range - s_step_pts[s_step_pts.size() - 2]);

    return new TSpline3(("spSRIM_" + particleKey).c_str(), s_step_pts.data(), y_step_pts.data(), s_step_pts.size(),
                        "b2,e2", 0, 0);
}

void convolutionSplineSRIMandGausDiffusion()
{
    // This macro is to check the convolution of the spline from SRIM and a gaussian diffusion
    // to see if it can reproduce the charge profile in the TPC

    auto srim = new ActPhysics::SRIM();
    srim->ReadTable("p", "../../Calibrations/SRIM/1H_900mb_CF4_95-5.txt");
    double range = 120.0;
    std::string particleKey = "p";
    double step = 0.5;

    TSpline3* sp = BuildSRIMspline(srim, range, particleKey, step);

    // Build histogram with the spline values
    int nbins = 1000;
    TH1D* hSpline = new TH1D("hSpline", "hSpline", nbins, 0, range);
    for(int i = 1; i <= nbins; i++)
    {
        double x = hSpline->GetBinCenter(i);
        double y = sp->Eval(x);
        hSpline->SetBinContent(i, y);
    }

    // copy the histogram to a new one to do the convolution
    TH1D* hConv = (TH1D*)hSpline->Clone("hConv");
    

    // Draw the spline histogram
    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
    hSpline->Draw();
}