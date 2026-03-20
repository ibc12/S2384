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

// ============================================================
// Convolve spline with Gaussian kernel to simulate diffusion (simple numerical convolution)
// ============================================================
TSpline3* BuildDiffusedSpline(TSpline3* spSRIM, double sigma, // sigma^2 = a + b*z
                                    int nPoints = 3000, double nSigma = 5.0)
{
    if(!spSRIM)
        return nullptr;

    double xmin = spSRIM->GetXmin();
    double xmax = spSRIM->GetXmax();

    double dx = (xmax - xmin) / (nPoints - 1);

    // --- sample the spline ---
    std::vector<double> x(nPoints), y(nPoints);
    for(int i = 0; i < nPoints; ++i)
    {
        x[i] = xmin + i * dx;
        y[i] = spSRIM->Eval(x[i]);
    }

    // --- extend range to allow for kernel tail ---
    int extraPoints = static_cast<int>(nSigma * sigma / dx);
    int nTotal = nPoints + extraPoints;

    std::vector<double> xExt(nTotal), yConv(nTotal, 0.0);

    for(int i = 0; i < nTotal; ++i)
        xExt[i] = xmin + i * dx;

    // --- causal smearing (only earlier points affect later positions) ---
    for(int i = 0; i < nTotal; ++i)
    {
        double sum = 0.0;
        double norm = 0.0;

        // only real SRIM points contribute
        for(int j = 0; j < nPoints; ++j)
        {
            if(x[j] > xExt[i])
                continue; // causality: only points <= current x contribute

            double d = xExt[i] - x[j];

            // limit kernel reach
            if(d > nSigma * sigma)
                continue;

            double w = std::exp(-0.5 * d * d / (sigma * sigma));

            sum += y[j] * w;
            norm += w;
        }

        if(norm > 0)
            yConv[i] = sum / norm;
    }

    // --- new spline ---
    auto spDiff = new TSpline3("spSRIM_diffused", xExt.data(), yConv.data(), nTotal);

    return spDiff;
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

TF1* FitSRIMtoChargeProfileFixedEnd(TH1* hCharge, TSpline3* spSRIM, const std::string& particleKey,
                                    double maxEndShift = 50.0) // cuanto puede moverse el final (mm)
{
    if(!spSRIM)
        return nullptr;

    // --- región útil del perfil ---
    double sOffset = FindPositionFromChargeFraction(hCharge, 0.02);
    double sEndData = FindPositionFromChargeFraction(hCharge, 0.98);
    double sEndFull = FindPositionFromChargeFraction(hCharge, 0.9999);

    std::cout << "Fit start offset = " << sOffset << " mm\n";
    std::cout << "Nominal end position = " << sEndFull << " mm\n";

    // --- final físico del spline (Bragg peak final) ---
    double rMax = spSRIM->GetXmax();

    std::string fname = "fSRIMfitFixedEnd_" + particleKey;

    TF1* f = new TF1(
        fname.c_str(),
        [spSRIM, sOffset, sEndFull, rMax](double* x, double* par)
        {
            double A = par[0];
            double deltaEnd = par[1];

            double s = x[0];

            if(s < sOffset)
                return 1e-9;

            // posición efectiva del final del track
            double sEndFit = sEndFull + deltaEnd;

            // alineación spline-datos
            double r = rMax - (sEndFit - s);

            if(r < spSRIM->GetXmin() || r > spSRIM->GetXmax())
                return 1e-9;

            return A * spSRIM->Eval(r);
        },
        sOffset, sEndData + maxEndShift, 2);

    // --- inicialización ---
    double integral = hCharge->Integral("width");
    double initAmp = integral / hCharge->GetNbinsX();

    f->SetParameters(initAmp, 0.0);

    f->SetParName(0, "Amplitude");
    f->SetParName(1, "EndShift");

    f->SetParLimits(0, 0, 1e12);
    f->SetParLimits(1, -maxEndShift, maxEndShift);

    f->SetNpx(3000);

    hCharge->Fit(f, "QR0", "", sOffset, sEndData + maxEndShift);

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
        model->SetBinError(ib, 1);
    }

    return model;
}

void fitTPCevent_diffusion() // LINEAR, CHI2 NORMALIZED, DQ/TL, RANGE
{
    PrettyStyle(false, false);
    // Get a charge profile histogram to do the fit
    // TFile* file = TFile::Open("./Outputs/hShifted_profile_lightD.root", "READ");
    TFile* file = TFile::Open("./Inputs/hProfile_experiment_IMPOSIBLE.root", "READ");
    // TFile* file = TFile::Open("./Inputs/hProfile_experiment_p_hard.root", "READ");
    // TFile* file = TFile::Open("../../Macros/L1PID/Outputs/qProfile_d_long.root", "READ"); // for sure d from
    // 7Li L1 PID
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

    // Check bin contents and errors for histogram
    for(int i = 1; i <= hShifted->GetNbinsX(); ++i)
    {
        double content = hShifted->GetBinContent(i);
        double error = hShifted->GetBinError(i);
        // std::cout << "Bin " << i << ": content = " << content << ", error = " << error << std::endl;
        hShifted->SetBinError(i, 1);
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

    auto* c = new TCanvas("cFit", "Fit TPC event", 800, 600);
    c->cd();
    hShifted->SetLineColor(kRed);
    hShifted->Draw("HISTE");


    std::map<std::string, TSpline3*> splineMap;
    for(const auto& key : particles)
    {
        splineMap[key] = BuildSRIMspline(
            srim, 250, key, 0.5); // If range of spline is bigger, the BP is badly reconstructed by the spline
    }

    double bestFitScore = 1e20;
    std::string bestFitParticle;

    for(const auto& key : particles)
    {
        std::cout << "----------------------------------------------\n";
        // auto fSpline = FitSRIMtoChargeProfileUniversal(hShifted, splineMap[key], key);
        auto fSpline = FitSRIMtoChargeProfileFixedEnd(hShifted, splineMap[key], key, 20);

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

        std::cout << "----------------------------------------------\n";
    }
    leg->Draw();
    std::cout << "================================================\n";
    std::cout << "Best fit WITHOUT diffusion: " << bestFitParticle << " (chi2=" << bestFitScore << ")\n";
    std::cout << "================================================\n";

    // Do the same for the diffused spline
    auto* cDiff = new TCanvas("cFitDiff", "Fit TPC event with diffusion", 800, 600);
    cDiff->cd();
    hShifted->SetLineColor(kRed);
    hShifted->Draw("HISTE");

    bestFitScore = 1e20;
    bestFitParticle = "";

    TLegend* leg1 = new TLegend(0.60, 0.65, 0.88, 0.88);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.035);

    for(const auto& key : particles)
    {
        auto spDiff = BuildDiffusedSpline(splineMap[key], 0.6, 2000, 5.0);
        auto fSplineDiff = FitSRIMtoChargeProfileFixedEnd(hShifted, spDiff, key + "_diff", 20);

        if(!fSplineDiff)
            continue;

        // ---------- style ----------
        fSplineDiff->SetLineColor(colorMap[key]);
        fSplineDiff->SetLineStyle(kDashed);
        fSplineDiff->SetLineWidth(3);
        fSplineDiff->Draw("same");

        leg1->AddEntry(fSplineDiff, (key + " diff fit").c_str(), "l");

        // ---------- ROOT fit chi2 ----------
        double chiFit = fSplineDiff->GetChisquare() / fSplineDiff->GetNDF();
        std::cout << "Particle " << key << " diff fit-chi2/NDF = " << chiFit << std::endl;
        if(chiFit < bestFitScore)
        {
            bestFitScore = chiFit;
            bestFitParticle = key;
        }

        std::cout << "----------------------------------------------\n";
    }
    leg1->Draw();
    std::cout << "================================================\n";
    std::cout << "Best fit WITH diffusion: " << bestFitParticle << " (chi2=" << bestFitScore << ")\n";
    std::cout << "================================================\n";


    // Detach histogram from file so it remains usable after closing the file/directory
    hShifted->SetDirectory(nullptr);
    file->Close();
}