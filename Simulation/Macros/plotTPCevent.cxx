#include "ActCluster.h"
#include "ActLine.h"
#include "ActSRIM.h"
#include "ActTPCParameters.h"
#include "ActVoxel.h"

#include <random>

#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TStyle.h>

#include <Math/Point3D.h>
#include <Math/Vector3D.h>
#include <cmath>
#include <map>
#include <tuple>
#include <vector>

using XYZPoint = ROOT::Math::XYZPoint;
using XYZVector = ROOT::Math::XYZVector;

// ============================================================
// Geometry
// ============================================================
constexpr double voxelSize = 1.0;           // mm
constexpr double Gmean = 3000.0;            // Mean gain
constexpr double theta = 0.7;                 // Polya parameter
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

    auto it = voxelMap.find(key);
    if(it == voxelMap.end())
    {
        ActRoot::Voxel v;
        v.SetPosition(ActRoot::Voxel::XYZPointF {float((ix + 0.5) * voxelSize), float((iy + 0.5) * voxelSize),
                                                 float((iz + 0.5) * voxelSize)});
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
                             std::map<voxelKey, ActRoot::Voxel>& voxelMap, std::vector<XYZPoint>& electrons)
{
    if(eLoss <= 0 || nPortions <= 0)
        return;

    const double W = 36.4 * 0.95 + 34.3 * 0.05; // eV / electron in 95% D2 + 5% CF4
    const double diffT = 0.10;                  // mm / sqrt(mm)
    const double diffL = 0;                     // mm / sqrt(mm) No diffusion in z direction

    double h = center.Y();
    if(h <= 0)
        return;

    double sigmaT = diffT * std::sqrt(h);
    double sigmaL = diffL * std::sqrt(h);

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
            electrons.emplace_back(x, y, z);
        }
    }
}

// ============================================================
// Divide track using SRIM
// ============================================================
void DivideTrackInSegments(ActPhysics::SRIM* srim, double range, const XYZVector& dirIn, const XYZPoint& rp,
                           double step, int nSub, std::map<voxelKey, ActRoot::Voxel>& voxelMap,
                           std::vector<XYZPoint>& electrons)
{
    XYZVector dir = dirIn.Unit();
    double E = srim->EvalInverse("light", range);

    for(double r = 0; r < range; r += step)
    {
        double Epost = srim->Slow("light", E, step);
        double eLoss = E - Epost;
        E = Epost;

        if(eLoss <= 0)
            continue;

        XYZPoint center = rp + dir * (r + 0.5 * step);
        if(center.Y() <= 0)
            continue;

        DivideSegmentInPortions(eLoss, nSub, center, voxelMap, electrons);
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
    cluster.SetVoxels(voxels);
    cluster.ReFit();

    ActRoot::Line line = cluster.GetLine();
    auto p0 = line.GetPoint();
    auto u = line.GetDirection().Unit();

    // --------------------------------------------------
    // Compute s-range
    // --------------------------------------------------
    double sMin = 1e9;
    double sMax = -1e9;

    for(const auto& v : voxels)
    {
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
    int nBins = std::max(1, int(std::ceil(length / ds)));

    // --------------------------------------------------
    // Histograms
    // --------------------------------------------------
    auto* prof =
        new TProfile("chargeProfileP", "Charge profile (mean);Track length [mm];Mean charge", nBins, sMin, sMax);

    auto* hist = new TH1D("chargeProfileH", "Charge profile (sum);Track length [mm];Charge", nBins, sMin, sMax);

    // --------------------------------------------------
    // Fill
    // --------------------------------------------------
    int nDiv = 3; // Divisions per axis
    double div = 1.0 / nDiv;

    for(const auto& v : voxels)
    {
        const auto& pos = v.GetPosition();
        double q = v.GetCharge();

        if(!subdivideVoxels)
        {
            // Standard: one point per voxel
            XYZVector d = XYZVector(pos) - XYZVector(p0);
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
                        XYZPoint miniPos(pos.X() + ix * div * voxelSize, pos.Y() + iy * div * voxelSize,
                                         pos.Z() + iz * div * voxelSize);

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
        hist->SetBinError(i, 0.0);

    hist->SetFillColorAlpha(kRed + 1, 0.35);
    hist->SetLineColor(kRed + 2);
    hist->SetLineWidth(2);

    return {prof, hist};
}


// ============================================================
// MAIN
// ============================================================
void plotTPCevent(double range = 120, double thetaDeg = 10, double phiDeg = -45)
{
    gStyle->SetOptStat(0);
    gRandom->SetSeed(0);

    ActRoot::TPCParameters tpc {"Actar"};

    auto* srim = new ActPhysics::SRIM;
    srim->ReadTable("light", "../../Calibrations/SRIM/1H_H2_CF4_95-5.txt");

    XYZPoint rp(tpc.X() / 2, tpc.Y() / 2, tpc.Z() / 2);

    double th = thetaDeg * TMath::DegToRad();
    double ph = phiDeg * TMath::DegToRad();
    XYZVector dir(std::cos(th), std::sin(th) * std::sin(ph), std::sin(th) * std::cos(ph));

    std::map<voxelKey, ActRoot::Voxel> voxelMap;
    std::vector<XYZPoint> electrons;

    DivideTrackInSegments(srim, range, dir, rp, 2.0, 5, voxelMap, electrons);

    // ================= Primary electrons plots =================
    TH2D* hXY = new TH2D("hXY", "XY;X [mm];Y [mm]", tpc.X(), 0, tpc.X(), tpc.Y(), 0, tpc.Y());

    TH2D* hXZ = new TH2D("hXZ", "XZ;X [mm];Z [mm]", tpc.X(), 0, tpc.X(), tpc.Z(), 0, tpc.Z());

    TH2D* hYZ = new TH2D("hYZ", "YZ;Z [mm];Y [mm]", tpc.Z(), 0, tpc.Z(), tpc.Y(), 0, tpc.Y());

    // ================= Charge plots =================
    TH2D* hXYq = new TH2D("hXYq", "Charge XY;X [mm];Y [mm]", tpc.X(), 0, tpc.X(), tpc.Y(), 0, tpc.Y());

    TH2D* hXZq = new TH2D("hXZq", "Charge XZ;X [mm];Z [mm]", tpc.X(), 0, tpc.X(), tpc.Z(), 0, tpc.Z());

    TH2D* hYZq = new TH2D("hYZq", "Charge YZ;Z [mm];Y [mm]", tpc.Z(), 0, tpc.Z(), tpc.Y(), 0, tpc.Y());

    // Fill the primary electrons plots
    for(const auto& e : electrons)
    {
        hXY->Fill(e.X(), e.Y());
        hXZ->Fill(e.X(), e.Z());
        hYZ->Fill(e.Z(), e.Y());
    }
    // Fill the charge plots
    for(const auto& [key, v] : voxelMap)
    {
        const auto& pos = v.GetPosition();
        double q = v.GetCharge();

        hXYq->Fill(pos.X(), pos.Y(), q);
        hXZq->Fill(pos.X(), pos.Z(), q);
        hYZq->Fill(pos.Z(), pos.Y(), q);
    }

    // ================= Profiles =================
    auto [profile, hist] = GetChargeProfile(voxelMap, true);

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
}
