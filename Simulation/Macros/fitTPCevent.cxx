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

// ============================================================
// Shift histogram along x-axis by a constant (for alignment) with x = 0 at the start of the track
// ============================================================
TH1* ShiftHistogram(TH1* h, double shift, const std::string& particleKey)
{
    int n = h->GetNbinsX();
    double xmin = h->GetBinLowEdge(1) + shift;
    double xmax = h->GetBinLowEdge(n + 1) + shift;

    TH1D* hnew = new TH1D(("shifted_" + particleKey).c_str(), h->GetTitle(), n, xmin, xmax);

    for(int i = 1; i <= n; i++)
        hnew->SetBinContent(i, h->GetBinContent(i));

    return hnew;
}

// ============================================================
// Integral with bin width (physically correct)
// ============================================================
double IntegralWidth(const TH1* h)
{
    return h->Integral("width");
}

// ============================================================
// Normalize histogram to integral = 1 (shape only)
// ============================================================
void NormalizeHistogram(TH1* h)
{
    double I = h->Integral("width");
    if(I > 0)
        h->Scale(1.0 / I);
}

// ============================================================
// Shape chi2 (robust for large charges)
// Does not use statistical errors -> compares shapes
// ============================================================
double Chi2Shape(const TH1* data, const TH1* model)
{
    double chi2 = 0.0;
    int n = 0;

    for(int i = 1; i <= data->GetNbinsX(); i++)
    {
        auto d = data->GetBinContent(i);
        auto m = model->GetBinContent(i);

        if(d <= 0 && m <= 0)
            continue;

        auto denom = d + m; // stable symmetric weight
        chi2 += (d - m) * (d - m) / denom;
        n++;
    }

    if(n == 0)
        return 1e12;

    return chi2 / n;
}

TSpline3* BuildSRIMspline(ActPhysics::SRIM* srim, double range, const std::string& particleKey, double step = 0.5,
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

    int nSteps = (int)s_step_pts.size();
    if(nSteps < 3)
    {
        std::cout << "Not enough points to build spline (nSteps=" << nSteps << "). Need at least 3 points.\n";
        return nullptr;
    }
    sp = new TSpline3(("spSRIM_" + particleKey).c_str(), s_step_pts.data(), y_step_pts.data(), nSteps, "b2,e2", 0, 0);

    return sp;
}

// ============================================================
// Find position along x-axis where cumulative integral reaches a given fraction of total charge
// ============================================================
double FindPositionFromChargeFraction(TH1* h, double frac)
{
    double total = h->Integral("width");
    double accum = 0.0;

    for(int i = 1; i <= h->GetNbinsX(); ++i)
    {
        accum += h->GetBinContent(i) * h->GetBinWidth(i);

        if(accum / total >= frac)
            return h->GetBinCenter(i);
    }

    return h->GetXaxis()->GetXmax();
}

TF1* FitSRIMtoChargeProfileUniversal(TH1* hCharge, TSpline3* spSRIM, const std::string& particleKey)
{
    if(!spSRIM)
        return nullptr;

    double sOffset = FindPositionFromChargeFraction(hCharge, 0.02);
    double sEnd = FindPositionFromChargeFraction(hCharge, 0.99);

    std::cout << "Fit start offset = " << sOffset << " mm\n";
    std::cout << "Fit end position = " << sEnd << " mm\n";

    double xmin = hCharge->GetXaxis()->GetXmin();
    double xmax = hCharge->GetXaxis()->GetXmax();

    std::string fname = "fSRIMfitUniversal_" + particleKey;

    TF1* f = new TF1(
        fname.c_str(),
        [spSRIM, sOffset, sEnd](double* x, double* par)
        {
            double A = par[0];
            double R = par[1];

            double s = x[0];

            if(s < sOffset)
                return 0.0;

            // Coordinate system change
            double r = R - (sEnd - s);

            if(r < spSRIM->GetXmin() || r > spSRIM->GetXmax())
                return 0.0;

            return A * spSRIM->Eval(r);
        },
        sOffset, sEnd, 2);

    double integral = hCharge->Integral("width");
    double initAmp = integral / hCharge->GetNbinsX();

    f->SetParameters(initAmp, 150.0);

    f->SetParName(0, "Amplitude");
    f->SetParName(1, "Range");

    f->SetParLimits(0, 0, 1e12);
    f->SetParLimits(1, 10, 400);

    f->SetNpx(500);

    hCharge->Fit(f, "RQ0", "", sOffset, sEnd);

    return f;
}

TH1* BuildModelHistogramFromTF1(const TH1* data, TF1* f, const std::string& name)
{
    TH1D* model = (TH1D*)data->Clone(name.c_str());
    model->Reset();

    for(int ib = 1; ib <= model->GetNbinsX(); ++ib)
    {
        double x = model->GetBinCenter(ib);
        double y = f->Eval(x);

        model->SetBinContent(ib, y);
        model->SetBinError(ib, std::sqrt(y));
    }

    return model;
}

void fitTPCevent()
{
    // Get a charge profile histogram to do the fit
    // TFile* file = TFile::Open("./Outputs/hShifted_profile_light.root", "READ");
    TFile* file = TFile::Open("./Inputs/hProfile_experiment.root", "READ");
    if(!file || file->IsZombie())
    {
        std::cerr << "Error: Could not open file ./Outputs/hShifted_profile_light.root" << std::endl;
        return;
    }

    // auto* hShifted = dynamic_cast<TH1D*>(file->Get("shifted_light"));
    auto* hShifted = dynamic_cast<TH1F*>(file->Get("hQProfile"));
    if(!hShifted)
    {
        std::cerr << "Error: Could not find histogram shifted_light in file ./Outputs/hShifted_profile_light.root"
                  << std::endl;
        std::cerr << "Contents of the file:\n";
        TIter next(file->GetListOfKeys());
        TKey* key;
        while((key = (TKey*)next()))
        {
            std::cout << " - " << key->GetName() << " (" << key->GetClassName() << ")\n";
        }
        file->Close();
        return;
    }

    // Do the fit here
    auto* srim = new ActPhysics::SRIM;
    srim->ReadTable("light", "../../Calibrations/SRIM/1H_900mb_CF4_95-5.txt");
    srim->ReadTable("lightD", "../../Calibrations/SRIM/2H_900mb_CF4_95-5.txt");
    srim->ReadTable("lightT", "../../Calibrations/SRIM/3H_900mb_CF4_95-5.txt");

    std::map<std::string, int> colorMap = {{"light", kBlue + 1}, {"lightD", kBlack}, {"lightT", kGreen + 2}};

    TLegend* leg = new TLegend(0.60, 0.65, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);

    std::vector<std::string> particles = {"light", "lightD", "lightT"};

    /////////////////////////////////////////////
    // Complex fit with range as free parameter
    /////////////////////////////////////////////


    auto* c = new TCanvas("cFit", "Fit TPC event", 800, 600);
    c->cd();
    hShifted->SetLineColor(kRed);
    hShifted->Draw("HIST");


    std::map<std::string, TSpline3*> splineMap;
    for(const auto& key : particles)
    {
        splineMap[key] = BuildSRIMspline(srim, 400, key, 0.5);
    }

    double bestScore = 1e12;
    double bestFitScore = 1e12;
    std::string bestParticle;
    std::string bestFitParticle;

    for(const auto& key : particles)
    {
        std::cout << "----------------------------------------------\n";
        auto fSpline = FitSRIMtoChargeProfileUniversal(hShifted, splineMap[key], key);

        if(!fSpline)
            continue;

        // ---------- style ----------
        fSpline->SetLineColor(colorMap[key]);
        fSpline->SetLineWidth(3);
        fSpline->Draw("same");

        leg->AddEntry(fSpline, (key + " fit").c_str(), "l");

        // ---------- ROOT fit chi2 ----------
        double chiFit = fSpline->GetChisquare() / fSpline->GetNDF();
        std::cout << "Particle " << key << "  fit-chi2/NDF = " << chiFit << std::endl;
        if(chiFit < bestFitScore)
        {
            bestFitScore = chiFit;
            bestFitParticle = key;
        }

        // ---------- build model histogram ----------
        TH1* model = BuildModelHistogramFromTF1(hShifted, fSpline, "model_" + key);

        // ---------- compare shape ----------
        TH1* dataNorm = (TH1*)hShifted->Clone(("dataNorm_" + key).c_str());
        TH1* modelNorm = (TH1*)model->Clone(("modelNorm_" + key).c_str());

        NormalizeHistogram(dataNorm);
        NormalizeHistogram(modelNorm);

        double chiShape = Chi2Shape(dataNorm, modelNorm);

        std::cout << "Particle " << key << "  shape-chi2 = " << chiShape << std::endl;

        if(chiShape < bestScore)
        {
            bestScore = chiShape;
            bestParticle = key;
        }
        std::cout << "----------------------------------------------\n";
    }
    leg->Draw();

    std::cout << "Best chi2 shape: " << bestParticle << " (chi2=" << bestScore << ")\n";
    std::cout << "Best usual chi2 fit: " << bestFitParticle << " (chi2=" << bestFitScore << ")\n";

    // Detach histogram from file so it remains usable after closing the file/directory
    hShifted->SetDirectory(nullptr);

    file->Close();
}