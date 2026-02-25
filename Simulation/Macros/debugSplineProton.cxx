#include "ActSRIM.h"

#include <TCanvas.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TMath.h>
#include <TSpline.h>
#include <TStyle.h>

#include <iostream>
#include <vector>

using namespace std;

// Standalone debug macro: build SRIM profile for "light" (proton) and
// draw both histogram-derived TGraph+TSpline and raw-sampled TGraph+TSpline
#include "ActSRIM.h"

#include <TCanvas.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TMath.h>
#include <TSpline.h>
#include <TStyle.h>

#include <iostream>
#include <vector>

using namespace std;

// Simplified debug macro for proton SRIM profile.
// - fills a histogram (hSRIM) with SRIM substeps (clamps negative E)
// - builds a histogram-based TGraph + TSpline
// - builds a sampled TGraph (no zero points) and a smooth monotone interpolant
// - prints concise diagnostics

void debugSplineProton(double range = 100.0, double step = 0.5, int nSubSteps = 10)
{
    gStyle->SetOptStat(0);

    // SRIM table
    ActPhysics::SRIM srim;
    srim.ReadTable("light", "../../Calibrations/SRIM/1H_900mb_CF4_95-5.txt");
    const string particle = "light";

    // Parameters
    const double sOffset = 5.0;
    const double ds = step / double(nSubSteps);

    // Histogram for SRIM (bin width ~= `step`)
    const int nBins = std::max(10, int(std::ceil((range + 10.0) / step)));
    const double sMin = 0.0;
    const double sMax = range + 10;
    TH1D hSRIM("dbg_hSRIM", "SRIM profile (hist);s [mm];dE [MeV]", nBins, sMin, sMax);

    // Integrate SRIM into histogram (clamp negative Epost)
    double E = srim.EvalInverse(particle, range);
    const double E_initial = E;
    bool ran_out = false;
    for(double r = 0; r < range; r += step)
    {
        double Epost = srim.Slow(particle, E, step);
        if(Epost < 0.0)
        {
            Epost = 0.0; // THIS WAS THE PROBLEM, I WAS NOT CHECKING FOR THIS
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
        const double dEsub = dE / nSubSteps;
        for(int i = 0; i < nSubSteps; ++i)
        {
            double s = r + (i + 0.5) * ds + sOffset;
            if(s >= sMin && s <= sMax)
                hSRIM.Fill(s, dEsub);
        }
        if(ran_out)
            break;
    }

    // Energy diagnostics
    const double E_after = E;
    const double sumContents = hSRIM.Integral();
    const double binW = hSRIM.GetBinWidth(1);
    const double sumWidthIntegral = hSRIM.Integral("width");
    const double diff = (E_initial - E_after) - sumContents;
    cout << "SRIM energy diagnostics:\n"
         << "  E_initial      = " << E_initial << " MeV\n"
         << "  E_after        = " << E_after << " MeV\n"
         << "  E_deposited    = " << (E_initial - E_after) << " MeV\n"
         << "  sum(hSRIM bins)= " << sumContents << " MeV (sum of bin contents)\n"
         << "  Integral(\"width\")= " << sumWidthIntegral << " (MeV*mm) [bin width=" << binW << " mm]\n"
         << "  difference     = " << diff << " MeV\n";

    // Histogram-based graph + spline (trim zeros)
    const int nb = hSRIM.GetNbinsX();
    vector<double> xb, yb;
    xb.reserve(nb);
    yb.reserve(nb);
    for(int i = 1; i <= nb; ++i)
    {
        xb.push_back(hSRIM.GetXaxis()->GetBinCenter(i));
        yb.push_back(hSRIM.GetBinContent(i));
    }
    // trim leading/trailing tiny bins
    const double tiny = 1e-12;
    int first = 0, last = int(xb.size()) - 1;
    while(first <= last && yb[first] <= tiny)
        ++first;
    while(last >= first && yb[last] <= tiny)
        --last;

    TGraph* gHist = nullptr;
    TSpline3* spHist = nullptr;
    if(last - first + 1 >= 1)
    {
        int m = last - first + 1;
        vector<double> xb_t(m), yb_t(m);
        for(int i = 0; i < m; ++i)
        {
            xb_t[i] = xb[first + i];
            yb_t[i] = yb[first + i];
        }
        gHist = new TGraph(m, xb_t.data(), yb_t.data());
        gHist->SetLineColor(kBlue + 1);
        gHist->SetLineWidth(2);
        gHist->SetMarkerStyle(20);
        gHist->SetMarkerSize(0.6);
        if(m >= 3)
        {
            spHist = new TSpline3("dbg_spHist", xb_t.data(), yb_t.data(), m, "b2,e2", 0, 0);
            spHist->SetLineColor(kBlue + 2);
            spHist->SetLineWidth(3);
        }
    }

    // Sampled profile: collect one point per integration step (user requested "only steps")
    vector<double> s_step_pts;
    s_step_pts.reserve(1024);
    vector<double> y_step_pts;
    y_step_pts.reserve(1024);

    E = srim.EvalInverse(particle, range);
    ran_out = false;
    for(double r = 0; r < range; r += step)
    {
        double Epost = srim.Slow(particle, E, step);
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
        double s = r + 0.5 * step + sOffset; // one point centered in the step
        s_step_pts.push_back(s);
        y_step_pts.push_back(dE);
        if(ran_out)
            break;
    }

    // Step diagnostics
    int stepN = (int)s_step_pts.size();
    double sumSteps = 0.0;
    double maxStep = 0.0;
    int imaxStep = 0;
    for(int i = 0; i < stepN; ++i)
    {
        sumSteps += y_step_pts[i];
        if(y_step_pts[i] > maxStep)
        {
            maxStep = y_step_pts[i];
            imaxStep = i;
        }
    }
    cout << "Sampled (per-step) profile diagnostics:\n  step points=" << stepN << "  sum(step deposits)=" << sumSteps
         << " MeV\n";
    if(stepN > 0)
        cout << "  max step point s=" << s_step_pts[imaxStep] << " value=" << maxStep << " MeV\n";

    // Build sampled TGraph from step points (exclude zero entries)
    TGraph* gSample = nullptr;
    if(!s_step_pts.empty())
    {
        gSample = new TGraph(stepN, s_step_pts.data(), y_step_pts.data());
        gSample->SetLineColor(kRed + 1);
        gSample->SetMarkerStyle(21);
    }

    // No smoothing or interpolation for step-sampled data per user request
    TGraph* gSampleSmooth = nullptr;

    // Draw
    TCanvas* c = new TCanvas("c_debugSpline", "Debug spline proton", 1000, 600);
    c->cd();
    hSRIM.SetLineColor(kBlue + 1);
    hSRIM.SetLineWidth(2);
    hSRIM.DrawClone("HIST");
    if(gHist)
        gHist->DrawClone("LP SAME");
    if(spHist)
        spHist->DrawClone("L SAME");
    double yMax = hSRIM.GetMaximum() * 4;
    hSRIM.GetYaxis()->SetRangeUser(0, yMax);
    if(gSample)
        gSample->DrawClone("LP SAME");
    if(gSampleSmooth)
        gSampleSmooth->DrawClone("L SAME");
    TLegend* leg = new TLegend(0.6, 0.65, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(&hSRIM, "hSRIM (hist)", "l");
    if(gHist)
        leg->AddEntry(gHist, "gHist (from hist bins)", "lp");
    if(spHist)
        leg->AddEntry(spHist, "spHist (spline from hist)", "l");
    if(gSample)
        leg->AddEntry(gSample, "gSample (raw sampled)", "lp");
    if(gSampleSmooth)
        leg->AddEntry(gSampleSmooth, "gSampleSmooth (monotone interp)", "l");
    leg->Draw();
}
