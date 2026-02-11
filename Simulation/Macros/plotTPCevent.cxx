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
constexpr double voxelSize = 2.0;               // mm
constexpr double Gmean = 3000.0;                // Mean gain
constexpr double theta = 0.7;                   // Polya parameter
constexpr double thresholdPadCharge = 5.4857e6; // that n electrons corresponds to 0.8789 pC
using voxelKey = std::tuple<int, int, int>;     // ix,iy,iz


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
                             std::map<voxelKey, ActRoot::Voxel>& voxelMap, std::vector<XYZPoint>& electrons)
{
    if(eLoss <= 0 || nPortions <= 0)
        return;

    const double W = 36.4 * 0.95 + 34.3 * 0.05; // eV / electron in 95% D2 + 5% CF4
    const double diffT = 0.06;                  // mm / sqrt(mm)- Aprox form data
    const double diffL = 0;                     // mm / sqrt(mm) No diffusion in z direction

    // Height in Z coordinate - distance to pad plane. Pad plane is in Z = 0 and asume rp is always in Z = 110. SO z
    // coordinate is distance to pad plane.
    double h = center.Z();
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
                           std::vector<XYZPoint>& electrons, bool isLight = true)
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
    auto* prof =
        new TProfile("chargeProfileP", "Charge profile (mean);Track length [mm];Mean charge", nBins, sMin - 5, sMax + 5);

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
                        XYZPoint miniPos(pos.X() * voxelSize + ix * div * voxelSize, pos.Y() * voxelSize + iy * div * voxelSize,
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
        hist->SetBinError(i, 0.0);

    hist->SetFillColorAlpha(kRed + 1, 0.35);
    hist->SetLineColor(kRed + 2);
    hist->SetLineWidth(2);

    return {prof, hist};
}


// ============================================================
// MAIN
// ============================================================
void plotTPCevent(double range = 120, double thetaDeg = 45, double phiDeg = -45)
{
    gStyle->SetOptStat(0);
    gRandom->SetSeed(0);

    ActRoot::TPCParameters tpc {"Actar"};

    auto* srim = new ActPhysics::SRIM;
    srim->ReadTable("light", "../../Calibrations/SRIM/1H_900mb_CF4_95-5.txt");
    srim->ReadTable("heavy", "../../Calibrations/SRIM/11Li_900mb_CF4_95-5.txt");

    XYZPoint rp(tpc.X() / 2, tpc.Y() / 2,
                110); // mm, starting point in the middle of the TPC for xy and 110 mm from the pad plane.

    double thLight = thetaDeg * TMath::DegToRad();
    double phLight = phiDeg * TMath::DegToRad();
    XYZVector dirLight(std::cos(thLight), std::sin(thLight) * std::sin(phLight), std::sin(thLight) * std::cos(phLight));

    double thHeavy = 0;
    double phHeavy = 0;
    XYZVector dirHeavy(std::cos(thHeavy), std::sin(thHeavy) * std::sin(phHeavy), std::sin(thHeavy) * std::cos(phHeavy));

    std::map<voxelKey, ActRoot::Voxel> voxelMapLight;
    std::map<voxelKey, ActRoot::Voxel> voxelMapHeavy;
    std::vector<XYZPoint> electronsLight;
    std::vector<XYZPoint> electronsHeavy;

    DivideTrackInSegments(srim, range, dirLight, rp, 2.0, 5, voxelMapLight, electronsLight, true);
    DivideTrackInSegments(srim, 3000, dirHeavy, rp, 2.0, 5, voxelMapHeavy, electronsHeavy, false);

    // ================= Primary electrons plots (units: mm) =================
    TH2D* hXY = new TH2D("hXY", "XY;X [mm];Y [mm]", tpc.X(), 0, tpc.X(), tpc.Y(), 0, tpc.Y());

    TH2D* hXZ = new TH2D("hXZ", "XZ;X [mm];Z [mm]", tpc.X(), 0, tpc.X(), tpc.Z(), 0, tpc.Z());

    TH2D* hYZ = new TH2D("hYZ", "YZ;Z [mm];Y [mm]", tpc.Z(), 0, tpc.Z(), tpc.Y(), 0, tpc.Y());
    // ================= Charge plots (units: voxels) =================
    TH2D* hXYq = new TH2D("hXYq", "Charge XY;X [vox];Y [vox]", tpc.X() / voxelSize, 0, tpc.X() / voxelSize,
                          tpc.Y() / voxelSize, 0, tpc.Y() / voxelSize);

    TH2D* hXZq = new TH2D("hXZq", "Charge XZ;X [vox];Z [vox]", tpc.X() / voxelSize, 0, tpc.X() / voxelSize,
                          tpc.Z() / voxelSize, 0, tpc.Z() / voxelSize);
    TH2D* hYZq = new TH2D("hYZq", "Charge YZ;Z [vox];Y [vox]", tpc.Z() / voxelSize, 0, tpc.Z() / voxelSize,
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

        if(q > 0) // here should be the threshold
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

        if(q > 0) // here should be the threshold
        {
            hXYq->Fill(pos.X(), pos.Y(), q);
            hXZq->Fill(pos.X(), pos.Z(), q);
            hYZq->Fill(pos.Z(), pos.Y(), q);
        }
    }

    // ================= Profiles =================
    auto [profile, hist] = GetChargeProfile(voxelMapLight, true);

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
