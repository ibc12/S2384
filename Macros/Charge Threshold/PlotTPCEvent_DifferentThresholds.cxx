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
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TSpline.h>
#include <TStyle.h>

#include <Math/Point3D.h>
#include <Math/Vector3D.h>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "../../PrettyStyle.C"

using XYZPoint = ROOT::Math::XYZPoint;
using XYZVector = ROOT::Math::XYZVector;

// ============================================================
// Geometry
// ============================================================
constexpr float voxelSize = 2.0;      // mm
ActRoot::TPCParameters tpc {"Actar"}; // TPC parameters
constexpr double Gmean = 3000.0;      // Mean gain
constexpr double theta = 0.7;         // Polya parameter
// constexpr double thresholdPadCharge = 5.4857e6; // that n electrons corresponds to 0.8789 pC
constexpr double thresholdPadCharge = 1e2; // that n electrons corresponds to 0.8789 pC
constexpr int yMinExclusionZone = 56;
constexpr int yMaxExclusionZone = 71;
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
                           std::vector<XYZPoint>& electrons, ActRoot::TPCParameters& tpc,
                           std::string particleType = "heavy")
{
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
std::pair<TGraph*, TH1D*> GetChargeProfile(const std::map<voxelKey, ActRoot::Voxel>& voxelMap, bool subdivideVoxels)
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
    auto* graph = new TGraph();
    graph->SetName("chargeProfileG");
    graph->SetTitle("Charge profile (mean);Track length [mm];Mean charge");

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
            v.SetPosition({(pos.X() + 0.5f) * voxelSize, (pos.Y() + 0.5f) * voxelSize,
                           (pos.Z() + 0.5f) * voxelSize}); // Convert voxel center from units of voxels to mm
            XYZVector d = XYZVector(v.GetPosition()) - XYZVector(p0);
            double s = d.Dot(u);

            // Get point in line with the same z coordinate than the voxel center
            // auto t {(v.GetPosition().Z() - p0.Z()) / u.Z()};
            // auto pointLineZcte = XYZPoint {p0.X() + u.X() * t, p0.Y() + u.Y() * t, p0.Z() + u.Z() * t};
            // Compute distance in s of this point
            // XYZVector d = XYZVector(v.GetPosition()) - XYZVector(p0);
            // double s = d.Dot(u);

            graph->SetPoint(graph->GetN(), s, q);
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

                        // auto t {(miniPos.Z() - p0.Z()) / u.Z()};
                        // auto pointLineZcte = XYZPoint {p0.X() + u.X() * t, p0.Y() + u.Y() * t, p0.Z() + u.Z() * t};
                        // XYZVector d = XYZVector(v.GetPosition()) - XYZVector(p0);
                        // double s = d.Dot(u);

                        graph->SetPoint(graph->GetN(), s, qSub);
                        hist->Fill(s, qSub);
                    }
                }
            }
        }
    }

    // for(int i = 1; i <= hist->GetNbinsX(); i++)
    //     hist->SetBinError(i, std::sqrt(hist->GetBinContent(i)));

    hist->SetFillColorAlpha(kRed + 1, 0.35);
    hist->SetLineColor(kRed + 2);
    hist->SetLineWidth(2);

    return {graph, hist};
}

// ============================================================
// Count pads outside exclusion zone (not taking into account angle or charge deposition)
// ============================================================
int PadsOutExclusionZone(const std::map<voxelKey, ActRoot::Voxel>& voxelMap1,
                         const std::map<voxelKey, ActRoot::Voxel>& voxelMap2, float threshold = 0)
{
    std::set<std::pair<int, int>> activePads;

    auto addPads = [&](const std::map<voxelKey, ActRoot::Voxel>& voxelMap)
    {
        for(const auto& [key, v] : voxelMap)
        {
            int ix = std::get<0>(key);
            int iy = std::get<1>(key);

            if((iy < yMinExclusionZone || iy > yMaxExclusionZone) && v.GetCharge() > threshold)
                activePads.insert({ix, iy});
        }
    };

    addPads(voxelMap1);
    addPads(voxelMap2);

    return activePads.size();
}

// ============================================================
// Shift histogram along x-axis by a constant (for alignment) with x = 0 at the start of the track
// ============================================================
TH1D* ShiftHistogram(TH1D* h, double shift, const std::string& particleKey)
{
    int n = h->GetNbinsX();
    double xmin = h->GetBinLowEdge(1) + shift;
    double xmax = h->GetBinLowEdge(n + 1) + shift;

    TH1D* hnew = new TH1D("hQProfile", h->GetTitle(), n, xmin, xmax);

    for(int i = 1; i <= n; i++)
    {
        hnew->SetBinContent(i, h->GetBinContent(i));
        hnew->SetBinError(i, h->GetBinError(i));
    }


    return hnew;
}

// ============================================================
// Result bundle for a single charge threshold: the three charge
// maps built with that threshold + the pads-out-of-exclusion-zone
// count. Building this does NOT require rerunning the SRIM /
// diffusion simulation, since it only re-filters the already
// simulated voxel maps.
// ============================================================
struct ThresholdResult
{
    double threshold = 0;
    int nPadsOutExclusionZone = 0;
    TH2D* hXYq = nullptr;
    TH2D* hXZq = nullptr;
    TH2D* hYZq = nullptr;
};

ThresholdResult
BuildThresholdPlots(double threshold, int index, const std::map<voxelKey, ActRoot::Voxel>& voxelMapLight,
                    const std::map<voxelKey, ActRoot::Voxel>& voxelMapHeavy, ActRoot::TPCParameters& tpc)
{
    ThresholdResult res;
    res.threshold = threshold;

    TString nameXY = TString::Format("hXYq_th%d", index);
    TString nameXZ = TString::Format("hXZq_th%d", index);
    TString nameYZ = TString::Format("hYZq_th%d", index);

    TString titleXY = TString::Format("Charge XY (thr = %.3g e^{-});X [pad];Y [pad]", threshold);
    TString titleXZ = TString::Format("Charge XZ (thr = %.3g e^{-});X [pad];Z [pad]", threshold);
    TString titleYZ = TString::Format("Charge YZ (thr = %.3g e^{-});Z [pad];Y [pad]", threshold);

    res.hXYq = new TH2D(nameXY, titleXY, tpc.X() / voxelSize, 0, tpc.X() / voxelSize, tpc.Y() / voxelSize, 0,
                        tpc.Y() / voxelSize);
    res.hXZq = new TH2D(nameXZ, titleXZ, tpc.X() / voxelSize, 0, tpc.X() / voxelSize, tpc.Z() / voxelSize, 0,
                        tpc.Z() / voxelSize);
    res.hYZq = new TH2D(nameYZ, titleYZ, tpc.Z() / voxelSize, 0, tpc.Z() / voxelSize, tpc.Y() / voxelSize, 0,
                        tpc.Y() / voxelSize);

    auto fillFrom = [&](const std::map<voxelKey, ActRoot::Voxel>& voxelMap)
    {
        for(const auto& [key, v] : voxelMap)
        {
            const auto& pos = v.GetPosition();
            double q = v.GetCharge();

            if(q > threshold) // here should be the threshold
            {
                res.hXYq->Fill(pos.X(), pos.Y(), q);
                res.hXZq->Fill(pos.X(), pos.Z(), q);
                res.hYZq->Fill(pos.Z(), pos.Y(), q);
            }
        }
    };
    fillFrom(voxelMapLight);
    fillFrom(voxelMapHeavy);

    res.nPadsOutExclusionZone = PadsOutExclusionZone(voxelMapLight, voxelMapHeavy, threshold);

    return res;
}

// ============================================================
// MAIN
// ============================================================
// chargeThresholds: list of charge thresholds (in electrons) to scan in a single
// run of the macro. For each value, the charge maps (XY/XZ/YZ) are plotted and the
// number of pads out of the exclusion zone is printed / summarized, without having
// to rerun the SRIM+diffusion simulation (which is the expensive part).
void PlotTPCEvent_DifferentThresholds()
{
    PrettyStyle();
    gRandom->SetSeed(0);

    // IMPORTANT: with several events processed in the same macro run, many histograms
    // built by the helper functions below share the same base name across events
    // (e.g. "chargeProfileH", "hXYq_th0", ...). By default ROOT auto-registers TH1-derived
    // histograms into the current directory and *deletes* any pre-existing object with the
    // same name, which would silently invalidate the previous event's canvases. Turning this
    // off makes every histogram a plain in-memory object, so per-event copies coexist safely.
    TH1::AddDirectory(kFALSE);

    auto* srim = new ActPhysics::SRIM;
    srim->ReadTable("light", "../../Calibrations/SRIM/1H_900mb_CF4_95-5.txt");
    srim->ReadTable("lightD", "../../Calibrations/SRIM/2H_900mb_CF4_95-5.txt");
    srim->ReadTable("lightT", "../../Calibrations/SRIM/3H_900mb_CF4_95-5.txt");
    srim->ReadTable("heavy", "../../Calibrations/SRIM/11Li_900mb_CF4_95-5.txt");

    // XYZPoint rp(tpc.X() / 2, tpc.Y() / 2,
    //             110); // mm, starting point in the middle of the TPC for xy and 110 mm from the pad plane.

    std::vector<double> chargeThresholds = {1e3, 1e4, 2e4, 3e4, 4e4, 5e4, 6e4, 7e4, 8e4, 9e4, 1e5}; // in electrons

    std::vector<double> thetaLights = {79, 77.2, 81.14, 149.5, 128.5, 149.4};
    std::vector<double> phiLights = {-122.8, 105.5, -100.6, -148.3, -87, -44};
    std::vector<double> thetaHeavy = {4.8, 5.5, 3.64, 2.8, 2.76, 1.36};
    std::vector<double> phiHeavy = {43.43, -72.73, 72.3, 30.57, 104.3, 96.3};
    std::vector<double> rp_x = {200.5, 16.4, 198., 208., 129.7, 149.8}; // in mm
    std::vector<double> rp_y = {128.5, 123.6, 126., 126., 123.6, 123.4}; // in mm
    std::vector<double> ranges = {30.8, 115.8, 23.4, 99.23, 136, 112.7};
    std::vector<std::string> lightStrings = {"lightD", "lightD", "lightD", "light", "light", "light"};

    // Event identifiers, same ordering/positions as the lists above (index i <-> eventNames[i])
    std::vector<std::string> eventNames = {"r66e292", "r66e973", "r66e1554", "r66e6583", "r66e8012", "r67e23132"};

    // Experimental pads-out-of-exclusion-zone count for each event, same ordering as eventNames.
    // Used below to build the simulated/experimental ratio vs threshold.
    std::vector<double> nPadsExperimental = {13, 278, 12, 50, 141, 100};

    // Collect the per-event "pads out of exclusion zone vs threshold" graphs so we can
    // also draw them together at the very end for a quick cross-event comparison.
    std::vector<TGraph*> gSummaryPerEvent;
    // Same, but for nPadsSimulated / nPadsExperimental vs threshold.
    std::vector<TGraph*> gRatioPerEvent;

    // ================= Loop over every event in the lists above =================
    for(int wantedIdx = 0; wantedIdx < (int)eventNames.size(); wantedIdx++)
    {
        const std::string& eventName = eventNames[wantedIdx];
        std::cout << "\n\n########## Processing event " << eventName << " (idx = " << wantedIdx << ") ##########"
                  << std::endl;

        std::string lightString = lightStrings[wantedIdx];

        double thLight = thetaLights[wantedIdx] * TMath::DegToRad();
        double phLight = phiLights[wantedIdx] * TMath::DegToRad();
        XYZVector dirLight(std::cos(thLight), std::sin(thLight) * std::sin(phLight),
                           std::sin(thLight) * std::cos(phLight));

        double range = ranges[wantedIdx];

        XYZPoint rp = {rp_x[wantedIdx], rp_y[wantedIdx], 110};

        double thHeavy = thetaHeavy[wantedIdx] * TMath::DegToRad();
        double phHeavy = phiHeavy[wantedIdx] * TMath::DegToRad();
        XYZVector dirHeavy(std::cos(thHeavy), std::sin(thHeavy) * std::sin(phHeavy),
                           std::sin(thHeavy) * std::cos(phHeavy));

        std::map<voxelKey, ActRoot::Voxel> voxelMapLight;
        std::map<voxelKey, ActRoot::Voxel> voxelMapHeavy;
        std::vector<XYZPoint> electronsLight;
        std::vector<XYZPoint> electronsHeavy;

        // ---- Expensive part: run ONCE per event regardless of how many thresholds are scanned ----
        DivideTrackInSegments(srim, range, dirLight, rp, 2, 5, voxelMapLight, electronsLight, tpc, lightString);
        DivideTrackInSegments(srim, 3000, dirHeavy, rp, 2, 5, voxelMapHeavy, electronsHeavy, tpc);

        // ================= Primary electrons plots (units: mm) =================
        TH2D* hXY = new TH2D(TString::Format("hXY_%s", eventName.c_str()),
                             TString::Format("XY (%s);X [mm];Y [mm]", eventName.c_str()), tpc.X(), 0, tpc.X(),
                             tpc.Y(), 0, tpc.Y());

        TH2D* hXZ = new TH2D(TString::Format("hXZ_%s", eventName.c_str()),
                             TString::Format("XZ (%s);X [mm];Z [mm]", eventName.c_str()), tpc.X(), 0, tpc.X(),
                             tpc.Z(), 0, tpc.Z());

        TH2D* hYZ = new TH2D(TString::Format("hYZ_%s", eventName.c_str()),
                             TString::Format("YZ (%s);Z [mm];Y [mm]", eventName.c_str()), tpc.Z(), 0, tpc.Z(),
                             tpc.Y(), 0, tpc.Y());

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

        // ================= Profiles (independent of the charge threshold) =================
        auto [graph, hist] = GetChargeProfile(voxelMapLight, true);
        graph->SetName(TString::Format("chargeProfileG_%s", eventName.c_str()));
        hist->SetName(TString::Format("chargeProfileH_%s", eventName.c_str()));

        // ================= Electrons canvas =================
        TCanvas* cEle = new TCanvas(TString::Format("cEle_%s", eventName.c_str()),
                                    TString::Format("Primary electrons - %s", eventName.c_str()), 1400, 450);
        cEle->Divide(3, 1);

        cEle->cd(1);
        hXY->Draw("COLZ");

        cEle->cd(2);
        hXZ->Draw("COLZ");

        cEle->cd(3);
        hYZ->Draw("COLZ");

        // ================= Charge profile canvas (independent of threshold) =================
        TCanvas* cProfile = new TCanvas(TString::Format("cProfile_%s", eventName.c_str()),
                                        TString::Format("Charge profile - %s", eventName.c_str()), 1500, 500);
        cProfile->Divide(3, 1);

        cProfile->cd(1);
        graph->SetMarkerStyle(4);
        graph->Draw("AP");

        cProfile->cd(2);
        hist->Draw("HIST");

        cProfile->cd(3);
        auto histProfileCopy = (TH1D*)hist->Clone(TString::Format("histProfileCopy_%s", eventName.c_str()));
        histProfileCopy->SetTitle(
            TString::Format("Charge profile shifted to origin (%s);Track length [mm];Charge", eventName.c_str()));
        double deltaS = -histProfileCopy->GetXaxis()->GetXmin();
        auto hShifted = ShiftHistogram(histProfileCopy, deltaS, "light");
        hShifted->SetName(TString::Format("hQProfile_%s", eventName.c_str()));
        hShifted->SetTitle(TString::Format(
            "Charge profile shifted to origin (%s);Track length from start of track [mm];Charge", eventName.c_str()));
        hShifted->Draw("HIST");

        // Save hShifted for later fit: create output file and write the histogram there.
        {
            TFile fout(
                TString::Format("./Outputs/hShifted_profile_%s_%s.root", eventName.c_str(), lightString.c_str()),
                "RECREATE");
            hShifted->Write();
            fout.Close();
        }
        // hShifted is owned by caller here; delete if you don't need it later.
        // delete hShifted;

        // ================= Scan over charge thresholds =================
        // For each threshold: build the charge maps (reusing the voxel maps already
        // simulated above), plot them, and report the pads out of exclusion zone.
        std::vector<ThresholdResult> results;
        results.reserve(chargeThresholds.size());

        for(size_t i = 0; i < chargeThresholds.size(); i++)
        {
            auto res = BuildThresholdPlots(chargeThresholds[i], (int)i, voxelMapLight, voxelMapHeavy, tpc);
            // BuildThresholdPlots only tags histograms by threshold index, so also stamp the
            // event name onto them for clarity when several events are processed in one run.
            res.hXYq->SetName(TString::Format("hXYq_th%zu_%s", i, eventName.c_str()));
            res.hXZq->SetName(TString::Format("hXZq_th%zu_%s", i, eventName.c_str()));
            res.hYZq->SetName(TString::Format("hYZq_th%zu_%s", i, eventName.c_str()));

            std::cout << "[" << eventName << "] Threshold = " << res.threshold
                      << " e-  ->  Pads out of exclusion zone = " << res.nPadsOutExclusionZone << std::endl;

            TString cName = TString::Format("cThr%zu_%s", i, eventName.c_str());
            TString cTitle = TString::Format("TPC event %s (thr = %.3g e-, pads out = %d)", eventName.c_str(),
                                             res.threshold, res.nPadsOutExclusionZone);

            TCanvas* cThr = new TCanvas(cName, cTitle, 1500, 500);
            cThr->Divide(3, 1);

            cThr->cd(1);
            res.hXYq->Draw("COLZ");
            cThr->cd(2);
            res.hXZq->Draw("COLZ");
            cThr->cd(3);
            res.hYZq->Draw("COLZ");

            results.push_back(res);
        }

        // ================= Summary: pads out of exclusion zone vs threshold (this event) =================
        std::cout << "\n===== Summary for " << eventName << ": pads out of exclusion zone vs threshold ====="
                  << std::endl;
        std::cout << std::left << std::setw(20) << "Threshold [e-]" << "Pads out of exclusion zone" << std::endl;
        for(const auto& res : results)
            std::cout << std::left << std::setw(20) << res.threshold << res.nPadsOutExclusionZone << std::endl;

        auto* gSummary = new TGraph();
        gSummary->SetName(TString::Format("gPadsVsThreshold_%s", eventName.c_str()));
        gSummary->SetTitle(TString::Format(
            "Pads out of exclusion zone vs threshold - %s;Charge threshold [e-];Pads out of exclusion zone",
            eventName.c_str()));
        for(const auto& res : results)
            gSummary->SetPoint(gSummary->GetN(), res.threshold, res.nPadsOutExclusionZone);

        TCanvas* cSummary = new TCanvas(TString::Format("cSummary_%s", eventName.c_str()),
                                        TString::Format("Pads vs threshold summary - %s", eventName.c_str()), 700,
                                        500);
        cSummary->SetLogx();
        gSummary->SetMarkerStyle(20);
        gSummary->SetMarkerColor(kAzure + 2);
        gSummary->SetLineColor(kAzure + 2);
        gSummary->Draw("APL");

        gSummaryPerEvent.push_back(gSummary);

        // ================= Simulated / experimental nPads ratio vs threshold (this event) =================
        double nExp = nPadsExperimental[wantedIdx];

        auto* gRatio = new TGraph();
        gRatio->SetName(TString::Format("gRatioVsThreshold_%s", eventName.c_str()));
        gRatio->SetTitle(TString::Format("Simu/Exp pads ratio vs threshold - %s (exp = %.0f);Charge threshold "
                                         "[e-];nPads_{simu} / nPads_{exp}",
                                         eventName.c_str(), nExp));
        for(const auto& res : results)
            gRatio->SetPoint(gRatio->GetN(), res.threshold, res.nPadsOutExclusionZone / nExp);

        TCanvas* cRatio = new TCanvas(TString::Format("cRatio_%s", eventName.c_str()),
                                      TString::Format("Simu/Exp pads ratio - %s", eventName.c_str()), 700, 500);
        cRatio->SetLogx();
        gRatio->SetMarkerStyle(20);
        gRatio->SetMarkerColor(kAzure + 2);
        gRatio->SetLineColor(kAzure + 2);
        gRatio->Draw("APL");

        gRatioPerEvent.push_back(gRatio);
    }

    // ================= Combined summary across all events =================
    TCanvas* cSummaryAll = new TCanvas("cSummaryAll", "Pads vs threshold summary - all events", 700, 500);
    cSummaryAll->SetLogx();
    auto* legAll = new TLegend(0.6, 0.6, 0.89, 0.89);
    std::vector<int> palette = {kAzure + 2, kRed + 1, kGreen + 2, kOrange + 1, kMagenta + 1, kBlack};

    // Drawing the first graph with "APL" and the rest with "PL SAME" only works if the
    // first graph's axes already cover every other graph's range: ROOT sets the axis
    // limits from whichever graph draws the frame, and later "SAME" graphs get clipped
    // to those limits instead of resizing them. So scan all graphs first and build an
    // explicit frame that spans the global x/y range before drawing any of them.
    double xMin = 1e300, xMax = -1e300, yMin = 1e300, yMax = -1e300;
    for(auto* g : gSummaryPerEvent)
    {
        double* xs = g->GetX();
        double* ys = g->GetY();
        for(int i = 0; i < g->GetN(); i++)
        {
            xMin = std::min(xMin, xs[i]);
            xMax = std::max(xMax, xs[i]);
            yMin = std::min(yMin, ys[i]);
            yMax = std::max(yMax, ys[i]);
        }
    }
    // Add a bit of headroom, especially on y (linear axis) so points don't sit on the edge.
    double yPad = 0.1 * (yMax - yMin);
    if(yPad <= 0)
        yPad = 1; // fallback if all events give the same nPads
    yMin = std::max(0.0, yMin - yPad);
    yMax = yMax + yPad;

    auto* frameAll = cSummaryAll->DrawFrame(xMin * 0.8, yMin, xMax * 1.2, yMax);
    frameAll->SetTitle("Pads out of exclusion zone vs threshold - all events;Charge threshold [e-];Pads out of "
                       "exclusion zone");

    for(size_t i = 0; i < gSummaryPerEvent.size(); i++)
    {
        auto* g = gSummaryPerEvent[i];
        g->SetMarkerStyle(20);
        g->SetMarkerColor(palette[i % palette.size()]);
        g->SetLineColor(palette[i % palette.size()]);
        g->Draw("PL SAME");
        legAll->AddEntry(g, eventNames[i].c_str(), "lp");
    }
    legAll->Draw();

    // ================= Combined simu/exp ratio summary across all events =================
    TCanvas* cRatioAll = new TCanvas("cRatioAll", "Simu/Exp pads ratio - all events", 700, 500);
    cRatioAll->SetLogx();
    auto* legRatioAll = new TLegend(0.6, 0.6, 0.89, 0.89);

    double rxMin = 1e300, rxMax = -1e300, ryMin = 1e300, ryMax = -1e300;
    for(auto* g : gRatioPerEvent)
    {
        double* xs = g->GetX();
        double* ys = g->GetY();
        for(int i = 0; i < g->GetN(); i++)
        {
            rxMin = std::min(rxMin, xs[i]);
            rxMax = std::max(rxMax, xs[i]);
            ryMin = std::min(ryMin, ys[i]);
            ryMax = std::max(ryMax, ys[i]);
        }
    }
    // Make sure the ratio = 1 reference line always fits inside the frame too.
    ryMin = std::min(ryMin, 1.0);
    ryMax = std::max(ryMax, 1.0);
    double ryPad = 0.1 * (ryMax - ryMin);
    if(ryPad <= 0)
        ryPad = 0.1;
    ryMin = std::max(0.0, ryMin - ryPad);
    ryMax = ryMax + ryPad;

    auto* frameRatioAll = cRatioAll->DrawFrame(rxMin * 0.8, ryMin, rxMax * 1.2, ryMax);
    frameRatioAll->SetTitle(
        "Simu/Exp pads ratio vs threshold - all events;Charge threshold [e-];nPads_{simu} / nPads_{exp}");

    for(size_t i = 0; i < gRatioPerEvent.size(); i++)
    {
        auto* g = gRatioPerEvent[i];
        g->SetMarkerStyle(20);
        g->SetMarkerColor(palette[i % palette.size()]);
        g->SetLineColor(palette[i % palette.size()]);
        g->Draw("PL SAME");
        legRatioAll->AddEntry(g, eventNames[i].c_str(), "lp");
    }

    // Reference line at ratio = 1 (perfect simu/exp agreement)
    auto* lineOne = new TLine(rxMin * 0.8, 1.0, rxMax * 1.2, 1.0);
    lineOne->SetLineStyle(2);
    lineOne->SetLineColor(kGray + 2);
    lineOne->Draw();

    legRatioAll->Draw();
}