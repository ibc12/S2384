#include "ActCluster.h"
#include "ActLine.h"
#include "ActSRIM.h"
#include "ActTPCParameters.h"
#include "ActVoxel.h"

#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TProfile.h>
#include <TRandom.h>
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
constexpr double voxelSize = 2.0;           // mm
using voxelKey = std::tuple<int, int, int>; // ix,iy,iz

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

    const double W = 30.0;     // eV / electron
    const double diffT = 0.10; // mm / sqrt(mm)
    const double diffL = 0.05; // mm / sqrt(mm)

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
            double y = gRandom->Gaus(center.Y(), sigmaL);
            double z = gRandom->Gaus(center.Z(), sigmaT);

            AddChargeToVoxel(x, y, z, 1.0, voxelMap);
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
std::pair<TProfile*, TH1D*> GetChargeProfile(const std::map<voxelKey, ActRoot::Voxel>& voxelMap)
{
    std::vector<ActRoot::Voxel> voxels;
    for(const auto& [k, v] : voxelMap)
        voxels.push_back(v);

    ActRoot::Cluster cluster;
    cluster.SetVoxels(voxels);
    cluster.ReFit();

    ActRoot::Line line = cluster.GetLine();
    auto p0 = line.GetPoint();
    auto u = line.GetDirection().Unit();

    double sMin = 1e9;
    double sMax = -1e9;

    for(const auto& v : voxels)
    {
        XYZVector d = XYZVector(v.GetPosition()) - XYZVector(p0);
        double s = d.Dot(u);
        sMin = std::min(sMin, s);
        sMax = std::max(sMax, s);
    }

    int nBins = 100;

    auto* prof =
        new TProfile("chargeProfileP", "Charge profile (mean);Track length [mm];Mean charge", nBins, sMin, sMax);

    auto* hist = new TH1D("chargeProfileH", "Charge profile (sum);Track length [mm];Charge", nBins, sMin, sMax);

    for(const auto& v : voxels)
    {
        XYZVector d = XYZVector(v.GetPosition()) - XYZVector(p0);
        double s = d.Dot(u);

        prof->Fill(s, v.GetCharge());
        hist->Fill(s, v.GetCharge());
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
void plotTPCevent(double range = 120, double thetaDeg = 45, double phiDeg = 45)
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

    // ================= Event plots =================
    TH2D* hXY = new TH2D("hXY", "XY;X [mm];Y [mm]", tpc.X(), 0, tpc.X(), tpc.Y(), 0, tpc.Y());

    TH2D* hXZ = new TH2D("hXZ", "XZ;X [mm];Z [mm]", tpc.X(), 0, tpc.X(), tpc.Z(), 0, tpc.Z());

    TH2D* hYZ = new TH2D("hYZ", "YZ;Z [mm];Y [mm]", tpc.Z(), 0, tpc.Z(), tpc.Y(), 0, tpc.Y());

    for(const auto& e : electrons)
    {
        hXY->Fill(e.X(), e.Y());
        hXZ->Fill(e.X(), e.Z());
        hYZ->Fill(e.Z(), e.Y());
    }

    // ================= Profiles =================
    auto [profile, hist] = GetChargeProfile(voxelMap);

    // ================= Canvas =================
    TCanvas* c = new TCanvas("c", "TPC event + charge profile", 1800, 800);
    c->Divide(3, 2);

    c->cd(1);
    hXY->Draw("COLZ");
    c->cd(2);
    hXZ->Draw("COLZ");
    c->cd(3);
    hYZ->Draw("COLZ");

    c->cd(4);
    profile->Draw();
    c->cd(5);
    hist->Draw("HIST");
}
