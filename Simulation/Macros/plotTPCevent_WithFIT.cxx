#include "ActCluster.h"
#include "ActLine.h"
#include "ActSRIM.h"
#include "ActTPCParameters.h"
#include "ActVoxel.h"

#include <random>

#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
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

using XYZPoint = ROOT::Math::XYZPoint;
using XYZVector = ROOT::Math::XYZVector;

// ============================================================
// Geometry
// ============================================================
constexpr double voxelSize = 2.0;               // mm
ActRoot::TPCParameters tpc {"Actar"};           // TPC parameters
constexpr double Gmean = 3000.0;                // Mean gain
constexpr double theta = 0.7;                   // Polya parameter
constexpr double thresholdPadCharge = 5.4857e6; // that n electrons corresponds to 0.8789 pC
constexpr int yMinExclusionZone = 55;
constexpr int yMaxExclusionZone = 70;
using voxelKey = std::tuple<int, int, int>; // ix,iy,iz


// ============================================================
// Polya function for gain distribution in Micromegas
// ============================================================
double Polya(double Gmean, double theta)
{
    static std::mt19937 gen(std::random_device {}());
    std::gamma_distribution<double> gamma(theta + 1.0, Gmean / (theta + 1.0));
    return gamma(gen);
}


// ============================================================
// Add charge to voxel map
// ============================================================
void AddChargeToVoxel(double x, double y, double z, double charge, std::map<voxelKey, ActRoot::Voxel>& voxelMap)
{
    int ix = std::floor(x / voxelSize);
    int iy = std::floor(y / voxelSize);
    int iz = std::floor(z / voxelSize);

    voxelKey key {ix, iy, iz};

    if(ix < 0 || ix >= (tpc.X() / voxelSize))
        return;
    if(iy < 0 || iy >= (tpc.Y() / voxelSize))
        return;

    auto it = voxelMap.find(key);
    if(it == voxelMap.end())
    {
        ActRoot::Voxel v;
        // Store voxel position in units of voxels as the corner of voxel
        v.SetPosition(ActRoot::Voxel::XYZPointF {float((ix)), float((iy)), float((iz))});
        v.SetCharge(charge);
        voxelMap.emplace(key, v);
    }
    else
    {
        it->second.SetCharge(it->second.GetCharge() + charge);
    }
}

// ============================================================
// Divide segment → electrons → voxels
// ============================================================
void DivideSegmentInPortions(double eLoss, int nPortions, const XYZPoint& center,
                             std::map<voxelKey, ActRoot::Voxel>& voxelMap, std::vector<XYZPoint>& electrons,
                             ActRoot::TPCParameters& tpc)
{
    if(eLoss <= 0 || nPortions <= 0)
        return;

    const double W = 36.4 * 0.95 + 34.3 * 0.05; // eV / electron in 95% D2 + 5% CF4
    const double diffT = 0.06;                  // mm / sqrt(mm) - Approx from data
    const double diffL = 0;                     // mm / sqrt(mm) No diffusion in z direction

    // Height in Z coordinate - distance to pad plane. Pad plane is at Z = 0 and rp is assumed to be at Z = 110.
    // Thus the z coordinate is the distance to the pad plane.
    double h = center.Z();
    if(h <= 0)
        return;

    double sigmaT = diffT * std::sqrt(h);
    double sigmaL = diffL * std::sqrt(h);

    // =====================================================
    // EARLY GEOMETRIC REJECTION (CRITICAL SPEEDUP)
    // =====================================================

    const double difXmax = tpc.X() + 20;
    const double difYmax = tpc.Y() + 20;

    if(center.X() < -20 || center.X() > difXmax || center.Y() < -20 || center.Y() > difYmax)
    {
        return; // entire cloud cannot reach pad plane
    }
    // =====================================================

    double portionE = eLoss / nPortions; // MeV
    double meanNe = portionE * 1e6 / W;

    for(int i = 0; i < nPortions; i++)
    {
        int nElectrons = gRandom->Poisson(meanNe);
        if(nElectrons <= 0)
            continue;


        for(int e = 0; e < nElectrons; e++)
        {
            double x = gRandom->Gaus(center.X(), sigmaT);
            double y = gRandom->Gaus(center.Y(), sigmaT);
            double z = gRandom->Gaus(center.Z(), sigmaL);

            double gain = Polya(Gmean, theta);
            AddChargeToVoxel(x, y, z, gain, voxelMap);
            // store electron positions in mm
            electrons.emplace_back(x, y, z);
        }
    }
}

// ============================================================
// Divide track using SRIM
// ============================================================
void DivideTrackInSegments(ActPhysics::SRIM* srim, double range, const XYZVector& dirIn, const XYZPoint& rp,
                           double step, int nSub, std::map<voxelKey, ActRoot::Voxel>& voxelMap,
                           std::vector<XYZPoint>& electrons, ActRoot::TPCParameters& tpc, bool isLight = true)
{
    std::string particleType = isLight ? "light" : "heavy";
    XYZVector dir = dirIn.Unit();
    double E = srim->EvalInverse(particleType, range);

    for(double r = 0; r < range; r += step)
    {
        double Epost = srim->Slow(particleType, E, step);
        double eLoss = E - Epost;
        E = Epost;

        if(eLoss <= 0)
            continue;

        XYZPoint center = rp + dir * (r + 0.5 * step);
        //  deltaZ > 0  -> track goes to larger Z (away from pad plane) - 146 mm to cathode
        //  deltaZ < 0  -> track goes to smaller Z (closer to pad plane) - 110 mm to pad plane
        // double DeltaZ = center.Z() - rp.Z();
        if(center.Z() <= 0 || center.Z() > 256)
            continue;

        DivideSegmentInPortions(eLoss, nSub, center, voxelMap, electrons, tpc);
    }
}

// ============================================================
// Build profile + histogram
// ============================================================
std::pair<TProfile*, TH1D*> GetChargeProfile(const std::map<voxelKey, ActRoot::Voxel>& voxelMap, bool subdivideVoxels)
{
    std::vector<ActRoot::Voxel> voxels;
    voxels.reserve(voxelMap.size());
    for(const auto& [k, v] : voxelMap)
        voxels.push_back(v);

    ActRoot::Cluster cluster;
    cluster.SetVoxels(voxels); // Voxels in pad units
    cluster.ReFit();

    ActRoot::Line line = cluster.GetLine();
    line.Scale(voxelSize, voxelSize); // Convert line parameters from voxel units to mm
    auto p0 = line.GetPoint();
    auto u = line.GetDirection().Unit();


    // --------------------------------------------------
    // Compute s-range
    // --------------------------------------------------
    double sMin = 1e9;
    double sMax = -1e9;

    for(auto v : voxels)
    {
        auto pos = v.GetPosition();
        v.SetPosition({pos.X() * voxelSize, pos.Y() * voxelSize,
                       pos.Z() * voxelSize}); // Convert voxel center from units of voxels to mm
        XYZVector d = XYZVector(v.GetPosition()) - XYZVector(p0);
        double s = d.Dot(u);
        sMin = std::min(sMin, s);
        sMax = std::max(sMax, s);
    }

    double length = sMax - sMin;
    if(length <= 0)
        length = voxelSize;

    // --------------------------------------------------
    // Adaptive bin size (voxel-limited resolution)
    // --------------------------------------------------
    double ds = voxelSize * std::max({std::abs(u.X()), std::abs(u.Y()), std::abs(u.Z())});
    std::cout << "Adaptive bin size: " << ds << " mm" << std::endl;

    ds = std::max(ds, 0.5 * voxelSize);
    int nBins = std::max(1, int(std::ceil((length + 10) / ds)));

    // --------------------------------------------------
    // Histograms
    // --------------------------------------------------
    auto* prof = new TProfile("chargeProfileP", "Charge profile (mean);Track length [mm];Mean charge", nBins, sMin - 5,
                              sMax + 5);

    auto* hist = new TH1D("chargeProfileH", "Charge profile (sum);Track length [mm];Charge", nBins, sMin - 5, sMax + 5);

    // --------------------------------------------------
    // Fill
    // --------------------------------------------------
    int nDiv = 3; // Divisions per axis
    double div = 1.0 / nDiv;

    for(auto v : voxels)
    {
        auto pos = v.GetPosition();
        double q = v.GetCharge();

        if(!subdivideVoxels)
        {
            // Standard: one point per voxel, take the center
            v.SetPosition({(pos.X() + 0.5) * voxelSize, (pos.Y() + 0.5) * voxelSize,
                           (pos.Z() + 0.5) * voxelSize}); // Convert voxel center from units of voxels to mm
            XYZVector d = XYZVector(v.GetPosition()) - XYZVector(p0);
            double s = d.Dot(u);

            prof->Fill(s, q);
            hist->Fill(s, q);
        }
        else
        {
            // Subdivide voxel into 3x3x3 mini-voxels
            double qSub = q / (nDiv * nDiv * nDiv);

            for(int ix = -1; ix <= 1; ix++)
            {
                for(int iy = -1; iy <= 1; iy++)
                {
                    for(int iz = -1; iz <= 1; iz++)
                    {
                        XYZPoint miniPos(pos.X() * voxelSize + ix * div * voxelSize,
                                         pos.Y() * voxelSize + iy * div * voxelSize,
                                         pos.Z() * voxelSize + iz * div * voxelSize);

                        XYZVector d = XYZVector(miniPos) - XYZVector(p0);
                        double s = d.Dot(u);

                        prof->Fill(s, qSub);
                        hist->Fill(s, qSub);
                    }
                }
            }
        }
    }

    for(int i = 1; i <= hist->GetNbinsX(); i++)
        hist->SetBinError(i, std::sqrt(hist->GetBinContent(i)));

    hist->SetFillColorAlpha(kRed + 1, 0.35);
    hist->SetLineColor(kRed + 2);
    hist->SetLineWidth(2);

    return {prof, hist};
}

// ============================================================
// Count pads outside exclusion zone (not taking into account angle or charge deposition)
// ============================================================
int PadsOutExclusionZone(const std::map<voxelKey, ActRoot::Voxel>& voxelMap1,
                         const std::map<voxelKey, ActRoot::Voxel>& voxelMap2)
{
    std::set<std::pair<int, int>> activePads;

    auto addPads = [&](const std::map<voxelKey, ActRoot::Voxel>& voxelMap)
    {
        for(const auto& [key, v] : voxelMap)
        {
            int ix = std::get<0>(key);
            int iy = std::get<1>(key);

            if(iy < yMinExclusionZone || iy > yMaxExclusionZone)
                activePads.insert({ix, iy});
        }
    };

    addPads(voxelMap1);
    addPads(voxelMap2);

    return activePads.size();
}

// ============================================================
// Integral with bin width (physically correct)
// ============================================================
double IntegralWidth(const TH1D* h)
{
    return h->Integral("width");
}

// ============================================================
// Normalize histogram to integral = 1 (shape only)
// ============================================================
void NormalizeHistogram(TH1D* h)
{
    double I = IntegralWidth(h);
    if(I > 0)
        h->Scale(1.0 / I);
}

// ============================================================
// Shape chi2 (robust for large charges)
// Does not use statistical errors -> compares shapes
// ============================================================
double Chi2Shape(const TH1D* data, const TH1D* model)
{
    double chi2 = 0.0;
    int n = 0;

    for(int i = 1; i <= data->GetNbinsX(); i++)
    {
        double d = data->GetBinContent(i);
        double m = model->GetBinContent(i);

        if(d <= 0 && m <= 0)
            continue;

        double denom = d + m; // stable symmetric weight
        chi2 += (d - m) * (d - m) / denom;
        n++;
    }

    if(n == 0)
        return 1e12;

    return chi2 / n;
}
// ============================================================
// Shift histogram along x-axis by a constant (for alignment) with x = 0 at the start of the track
// ============================================================
TH1D* ShiftHistogram(TH1D* h, double shift, const std::string& particleKey)
{
    int n = h->GetNbinsX();
    double xmin = h->GetBinLowEdge(1) + shift;
    double xmax = h->GetBinLowEdge(n + 1) + shift;

    TH1D* hnew = new TH1D(("shifted_" + particleKey).c_str(), h->GetTitle(), n, xmin, xmax);

    for(int i = 1; i <= n; i++)
        hnew->SetBinContent(i, h->GetBinContent(i));

    return hnew;
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
// Find offset along x-axis where cumulative integral reaches a given fraction of total charge (for alignment)
// ============================================================
double FindOffsetFromChargeFraction(const TH1D* h, double frac = 0.02)
{
    double total = h->Integral("width");
    if(total <= 0)
        return h->GetXaxis()->GetXmin();

    double accum = 0.0;

    for(int i = 1; i <= h->GetNbinsX(); ++i)
    {
        accum += h->GetBinContent(i) * h->GetBinWidth(i);

        if(accum >= frac * total)
            return h->GetBinLowEdge(i);
    }

    return h->GetXaxis()->GetXmin();
}

TF1* FitSRIMtoChargeProfile(ActPhysics::SRIM* srim, double range, TH1D* hCharge, const std::string& particleKey,
                            double step = 0.5)
{
    // --- align charge axis to start at 0
    double deltaS = -hCharge->GetXaxis()->GetXmin();
    // Moving the axis outside the function call does not change the histogram bin values
    // So we have to shift manually the histogram values
    auto hShifted = ShiftHistogram(hCharge, deltaS, particleKey);

    // compute spline
    double sOffset = 5.0;

    // ---------------------------------------------------------
    // Build a per-step sampled TGraph (one point per full integration step)
    // ---------------------------------------------------------

    // Build TGraph and spline from per-step points
    std::string spname = "spSRIM_" + particleKey;
    TSpline3* spSRIM = nullptr;
    spSRIM = BuildSRIMspline(srim, range, particleKey, step, sOffset);


    TF1* fSpline = nullptr;
    if(spSRIM)
    {
        std::string fname = "fSRIMfit_" + particleKey;

        double fitMin = spSRIM->GetXmin();
        double fitMax = spSRIM->GetXmax();

        fSpline = new TF1(
            fname.c_str(),
            [spSRIM, fitMin, fitMax](double* x, double* par)
            {
                if(x[0] < fitMin || x[0] > fitMax)
                    return 0.0;
                return par[0] * spSRIM->Eval(x[0]);
            },
            fitMin, fitMax, 1);

        // Initial amplitude is related with integrals of hist and spline
        double histIntegral = hCharge->Integral("width");
        double initAmp = histIntegral / hCharge->GetNbinsX();
        // double splineIntegral = 0;
        // double xMin = spSRIM->GetXmin();
        // double xMax = spSRIM->GetXmax();
        // double dx = (xMax - xMin) / (range / step - 1);
        // for(int i = 0; i < range / step - 1; i++)
        // {
        //     double y1 = spSRIM->Eval(xMin + i * dx);
        //     double y2 = spSRIM->Eval(xMin + (i + 1) * dx);
        //     splineIntegral += 0.5 * (y1 + y2) * dx;
        // }
        // double initAmp = histIntegral / splineIntegral;

        fSpline->SetNpx(400);
        fSpline->SetParameter(0, initAmp);
        fSpline->SetParLimits(0, 0, 1e12);
        fSpline->SetParameter(1, 5);

        hShifted->Fit(fSpline, "MRV0", "", fitMin, fitMax);

        double amp = fSpline->GetParameter(0);
        double ampErr = fSpline->GetParError(0);
        std::cout << "SRIM spline fit amplitude (" << particleKey << ") = " << amp << " +/- " << ampErr
                  << " electrons / MeV\n";
    }
    else
    {
        std::cout << "Not enough per-step points for spline fit (nSteps=" << range / step
                  << "). Skipping spline fit.\n";
    }

    return fSpline;
}

TH1D* BuildModelHistogramFromTF1(const TH1D* data, TF1* f, const std::string& name)
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


// ============================================================
// MAIN
// ============================================================
void plotTPCevent_WithFIT(double range = 120, double thetaDeg = 45, double phiDeg = -45)
{
    gStyle->SetOptStat(0);
    gRandom->SetSeed(0);

    // double chargeThreshold = thresholdPadCharge; // Threshold in electrons
    double chargeThreshold = 0;

    auto* srim = new ActPhysics::SRIM;
    srim->ReadTable("light", "../../Calibrations/SRIM/1H_900mb_CF4_95-5.txt");
    srim->ReadTable("lightD", "../../Calibrations/SRIM/2H_900mb_CF4_95-5.txt");
    srim->ReadTable("lightT", "../../Calibrations/SRIM/3H_900mb_CF4_95-5.txt");
    srim->ReadTable("heavy", "../../Calibrations/SRIM/11Li_900mb_CF4_95-5.txt");

    XYZPoint rp(tpc.X() / 2, tpc.Y() / 2,
                110); // mm, starting point in the middle of the TPC for xy and 110 mm from the pad plane.

    double thLight = thetaDeg * TMath::DegToRad();
    double phLight = phiDeg * TMath::DegToRad();
    XYZVector dirLight(std::cos(thLight), std::sin(thLight) * std::sin(phLight), std::sin(thLight) * std::cos(phLight));

    double thHeavy = 5 * TMath::DegToRad();
    double phHeavy = 10 * TMath::DegToRad();
    XYZVector dirHeavy(std::cos(thHeavy), std::sin(thHeavy) * std::sin(phHeavy), std::sin(thHeavy) * std::cos(phHeavy));

    std::map<voxelKey, ActRoot::Voxel> voxelMapLight;
    std::map<voxelKey, ActRoot::Voxel> voxelMapHeavy;
    std::vector<XYZPoint> electronsLight;
    std::vector<XYZPoint> electronsHeavy;

    DivideTrackInSegments(srim, range, dirLight, rp, 2.0, 5, voxelMapLight, electronsLight, tpc, true);
    DivideTrackInSegments(srim, 3000, dirHeavy, rp, 2.0, 5, voxelMapHeavy, electronsHeavy, tpc, false);

    // ================= Primary electrons plots (units: mm) =================
    TH2D* hXY = new TH2D("hXY", "XY;X [mm];Y [mm]", tpc.X(), 0, tpc.X(), tpc.Y(), 0, tpc.Y());

    TH2D* hXZ = new TH2D("hXZ", "XZ;X [mm];Z [mm]", tpc.X(), 0, tpc.X(), tpc.Z(), 0, tpc.Z());

    TH2D* hYZ = new TH2D("hYZ", "YZ;Z [mm];Y [mm]", tpc.Z(), 0, tpc.Z(), tpc.Y(), 0, tpc.Y());
    // ================= Charge plots (units: voxels) =================
    TH2D* hXYq = new TH2D("hXYq", "Charge XY;X [pad];Y [pad]", tpc.X() / voxelSize, 0, tpc.X() / voxelSize,
                          tpc.Y() / voxelSize, 0, tpc.Y() / voxelSize);

    TH2D* hXZq = new TH2D("hXZq", "Charge XZ;X [pad];Z [pad]", tpc.X() / voxelSize, 0, tpc.X() / voxelSize,
                          tpc.Z() / voxelSize, 0, tpc.Z() / voxelSize);
    TH2D* hYZq = new TH2D("hYZq", "Charge YZ;Z [pad];Y [pad]", tpc.Z() / voxelSize, 0, tpc.Z() / voxelSize,
                          tpc.Y() / voxelSize, 0, tpc.Y() / voxelSize);
    // Fill the primary electrons plots
    for(const auto& e : electronsLight)
    {
        hXY->Fill(e.X(), e.Y());
        hXZ->Fill(e.X(), e.Z());
        hYZ->Fill(e.Z(), e.Y());
    }
    for(const auto& e : electronsHeavy)
    {
        hXY->Fill(e.X(), e.Y());
        hXZ->Fill(e.X(), e.Z());
        hYZ->Fill(e.Z(), e.Y());
    }
    // Fill the charge plots
    for(const auto& [key, v] : voxelMapLight)
    {
        const auto& pos = v.GetPosition();
        double q = v.GetCharge();

        if(q > chargeThreshold) // here should be the threshold
        {
            hXYq->Fill(pos.X(), pos.Y(), q);
            hXZq->Fill(pos.X(), pos.Z(), q);
            hYZq->Fill(pos.Z(), pos.Y(), q);
        }
    }
    for(const auto& [key, v] : voxelMapHeavy)
    {
        const auto& pos = v.GetPosition();
        double q = v.GetCharge();

        if(q > chargeThreshold) // here should be the threshold
        {
            hXYq->Fill(pos.X(), pos.Y(), q);
            hXZq->Fill(pos.X(), pos.Z(), q);
            hYZq->Fill(pos.Z(), pos.Y(), q);
        }
    }

    // ================= Profiles =================
    auto [profile, hist] = GetChargeProfile(voxelMapLight, true);

    // =============== Pads out of exclusion zone =================
    int nPadsOutExclusionZone = PadsOutExclusionZone(voxelMapLight, voxelMapHeavy);
    std::cout << "Number of pads out of exclusion zone: " << nPadsOutExclusionZone << std::endl;

    // ================= Electrons canvas =================
    TCanvas* cEle = new TCanvas("cEle", "Primary electrons", 1400, 450);
    cEle->Divide(3, 1);

    cEle->cd(1);
    hXY->Draw("COLZ");

    cEle->cd(2);
    hXZ->Draw("COLZ");

    cEle->cd(3);
    hYZ->Draw("COLZ");

    // ================= Main canvas =================
    TCanvas* c = new TCanvas("c", "TPC event + charge profile", 1800, 800);
    c->Divide(3, 2);

    c->cd(1);
    hXYq->Draw("COLZ");
    c->cd(2);
    hXZq->Draw("COLZ");
    c->cd(3);
    hYZq->Draw("COLZ");
    c->cd(4);
    profile->Draw();
    c->cd(5);
    hist->Draw("HIST");
    c->cd(6);
    auto histProfileCopy = (TH1D*)hist->Clone("histProfileCopy");
    double deltaS = -histProfileCopy->GetXaxis()->GetXmin();
    histProfileCopy->GetXaxis()->SetLimits(histProfileCopy->GetXaxis()->GetXmin() + deltaS,
                                           histProfileCopy->GetXaxis()->GetXmax() + deltaS);
    histProfileCopy->SetTitle("Charge profile with SRIM fit;Track length [mm];Charge");
    histProfileCopy->Draw("HIST");
    // FitSRIMtoChargeProfile(srim, range, histProfileCopy, "light", 0.3);

    std::map<std::string, int> colorMap = {{"light", kBlue + 1}, {"lightD", kBlack}, {"lightT", kGreen + 2}};

    TLegend* leg = new TLegend(0.60, 0.65, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);

    std::vector<std::string> particles = {"light", "lightD", "lightT"};

    double bestScore = 1e12;
    std::string bestParticle;

    double bestFitScore = 1e12;
    std::string bestFitParticle;

    for(const auto& key : particles)
    {
        auto fSpline = FitSRIMtoChargeProfile(srim, range, histProfileCopy, key, 0.5);

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
        TH1D* model = BuildModelHistogramFromTF1(histProfileCopy, fSpline, "model_" + key);

        // ---------- compare shape ----------
        TH1D* dataNorm = (TH1D*)histProfileCopy->Clone(("dataNorm_" + key).c_str());
        TH1D* modelNorm = (TH1D*)model->Clone(("modelNorm_" + key).c_str());

        NormalizeHistogram(dataNorm);
        NormalizeHistogram(modelNorm);

        double chiShape = Chi2Shape(dataNorm, modelNorm);

        std::cout << "Particle " << key << "  shape-chi2 = " << chiShape << std::endl;

        if(chiShape < bestScore)
        {
            bestScore = chiShape;
            bestParticle = key;
        }
    }
    leg->Draw();

    std::cout << "Best chi2 shape: " << bestParticle << " (chi2=" << bestScore << ")\n";
    std::cout << "Best usual chi2 fit: " << bestFitParticle << " (chi2=" << bestFitScore << ")\n";
}
