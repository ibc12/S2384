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
double Chi2Shape(const TH1* data, const TH1* model, const TF1* f)
{
    double xmin = f->GetXmin();
    double xmax = f->GetXmax();

    int binMin = data->GetXaxis()->FindBin(xmin);
    int binMax = data->GetXaxis()->FindBin(xmax);

    // --- integrales SOLO en rango del fit ---
    double intData = 0.0;
    double intModel = 0.0;

    for(int i = binMin; i <= binMax; i++)
    {
        double bw = data->GetBinWidth(i);
        intData += data->GetBinContent(i) * bw;
        intModel += model->GetBinContent(i) * bw;
    }

    if(intData <= 0 || intModel <= 0)
        return 1e12;

    // --- chi2 de forma ---
    double chi2 = 0.0;
    int n = 0;

    for(int i = binMin; i <= binMax; i++)
    {
        double d = data->GetBinContent(i) / intData;
        double m = model->GetBinContent(i) / intModel;

        if(d <= 0 && m <= 0)
            continue;

        double denom = d + m;
        chi2 += (d - m) * (d - m) / denom;
        n++;
    }

    if(n == 0)
        return 1e12;

    return chi2 / n;
}

double Chi2Manually(const TH1* data, const TH1* model, const TF1* f)
{
    double xmin = f->GetXmin();
    double xmax = f->GetXmax();

    int binMin = data->GetXaxis()->FindBin(xmin);
    int binMax = data->GetXaxis()->FindBin(xmax);

    double chi2 = 0.0;
    int n = 0;

    for(int i = binMin; i <= binMax; i++)
    {
        double d = data->GetBinContent(i);
        double m = model->GetBinContent(i);

        double ed = data->GetBinError(i);
        double em = model->GetBinError(i);

        // double sigma2 = ed * ed + em * em;
        double sigma2 = 1.0;

        if(sigma2 <= 0)
            continue;

        chi2 += (d - m) * (d - m) / sigma2;
        n++;
    }

    if(n == 0)
        return 1e12;

    return chi2 / n; // reduced chi2
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

    // std::cout << "Spline max: " << sp->GetXmax() << "\n";
    // std::cout << "Value near max: " << sp->Eval(range - 1) << "\n";

    return sp;
}

TGraph* BuildSRIMgraph(ActPhysics::SRIM* srim, double range, const std::string& particleKey, double step = 0.6,
                       double sOffset = 0.0)
{
    TGraph* gr = nullptr;
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

    if(!y_step_pts.empty() && y_step_pts.back() > 0.0)
    {
        s_step_pts.push_back(range);
        y_step_pts.push_back(0.0);
    }

    int nSteps = (int)s_step_pts.size();
    if(nSteps < 2) // para lineal solo necesitas 2
    {
        std::cout << "Not enough points to build graph (nSteps=" << nSteps << "). Need at least 2 points.\n";
        return nullptr;
    }

    gr = new TGraph(nSteps, s_step_pts.data(), y_step_pts.data());
    gr->SetName(("grSRIM_" + particleKey).c_str());

    return gr;
}

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
                                    double maxEndShift = 20.0) // cuanto puede moverse el final (mm)
{
    if(!spSRIM)
        return nullptr;

    // --- región útil del perfil ---
    double sOffset = FindPositionFromChargeFraction(hCharge, 0.02);
    double sEndData = FindPositionFromChargeFraction(hCharge, 0.98);
    double sEndFull = FindPositionFromChargeFraction(hCharge, 0.9999);

    // std::cout << "Fit start offset = " << sOffset << " mm\n";
    // std::cout << "Nominal end position = " << sEndData << " mm\n";

    // --- final físico del spline (Bragg peak final) ---
    double rMax = spSRIM->GetXmax();

    std::string fname = "fSRIMfitFixedEnd_" + particleKey;

    TF1* f = new TF1(
        fname.c_str(),
        [spSRIM, sOffset, sEndData, rMax](double* x, double* par)
        {
            double A = par[0];
            double deltaEnd = par[1];

            double s = x[0];

            if(s < sOffset)
                return 1e-9;

            // posición efectiva del final del track
            double sEndFit = sEndData + deltaEnd;

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

TF1* FitSRIMtoChargeProfileFixedEndGraph(TH1* hCharge, TGraph* grSRIM, const std::string& particleKey,
                                         double maxEndShift = 20.0) // cuanto puede moverse el final (mm)
{
    if(!grSRIM)
        return nullptr;

    // --- región útil del perfil ---
    double sOffset = FindPositionFromChargeFraction(hCharge, 0.02);
    double sEndData = FindPositionFromChargeFraction(hCharge, 0.98);
    double sEndFull = FindPositionFromChargeFraction(hCharge, 0.9999);

    // std::cout << "Fit start offset = " << sOffset << " mm\n";
    // std::cout << "Nominal end position = " << sEndData << " mm\n";

    // --- final físico del spline (Bragg peak final) ---
    double rMax = grSRIM->GetX()[grSRIM->GetN() - 1];

    std::string fname = "fSRIMfitFixedEnd_" + particleKey;

    TF1* f = new TF1(
        fname.c_str(),
        [grSRIM, sOffset, sEndData, rMax](double* x, double* par)
        {
            double A = par[0];
            double deltaEnd = par[1];

            double s = x[0];

            if(s < sOffset)
                return 1e-9;

            // posición efectiva del final del track
            double sEndFit = sEndData + deltaEnd;

            // alineación spline-datos
            double r = rMax - (sEndFit - s);

            if(r < grSRIM->GetX()[0] || r > grSRIM->GetX()[grSRIM->GetN() - 1])
                return 1e-9;

            return A * grSRIM->Eval(r);
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

void checkGoodFitsBraggPeak()
{
    PrettyStyle(false, true);

    // Get all L1 experimental data for 7Li; cut in L1 and filter bad events, select elastic, check good fits to
    // deuterium

    // Get all data for 7Li (there is no triton there, so easier to see deuterium)
    std::string beam {"7Li"};
    std::string target {"d"};
    std::string light {"p"};

    std::string dataconf {};
    if(beam == "11Li")
        dataconf = "../../configs/data_11Li.conf";
    else if(beam == "7Li")
        dataconf = "../../configs/data_7Li.conf";
    else
        throw std::runtime_error("Beam cannot differ from 11Li or 7Li");

    // Read data
    ActRoot::DataManager dataman {dataconf, ActRoot::ModeType::EMerge};
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    auto chain3 {dataman.GetChain(ActRoot::ModeType::EFilter)};
    chain->AddFriend(chain2.get());
    chain->AddFriend(chain3.get());

    // Get drift parameter
    ActRoot::InputParser parser {};
    parser.ReadFile("../../configs/detector.conf");
    auto driftBlock = parser.GetBlock("Merger");
    auto driftFactor = driftBlock->GetDouble("DriftFactor"); // in mm^2/us

    // Get SRIM files to do the fits
    auto* srim = new ActPhysics::SRIM;
    srim->ReadTable("light", "../../Calibrations/SRIM/1H_900mb_CF4_95-5.txt");
    srim->ReadTable("lightD", "../../Calibrations/SRIM/2H_900mb_CF4_95-5.txt");
    srim->ReadTable("lightT", "../../Calibrations/SRIM/3H_900mb_CF4_95-5.txt");
    srim->ReadTable("light3He", "../../Calibrations/SRIM/3He_900mb_CF4_95-5.txt");
    srim->ReadTable("light4He", "../../Calibrations/SRIM/4He_900mb_CF4_95-5.txt");

    // Create the splines
    std::vector<std::string> particles = {"light", "lightD", "lightT"};
    std::map<std::string, TSpline3*> splineMap;
    for(const auto& key : particles)
    {
        splineMap[key] = BuildSRIMspline(
            srim, 250, key, 0.5); // If range of spline is bigger, the BP is badly reconstructed by the spline
    }
    std::map<std::string, TGraph*> graphMap;
    for(const auto& key : particles)
    {
        graphMap[key] = BuildSRIMgraph(srim, 250, key, 0.5);
    }

    // RDataFrame
    // ROOT::EnableImplicitMT();
    ROOT::RDataFrame dforigin {*chain};

    // Filter GATCONF == L1
    auto df =
        dforigin.Filter([](ActRoot::ModularData& mod, ActRoot::MergerData& mer)
                        { return mod.Get("GATCONF") == 8 && (mer.fLightIdx != -1); }, {"ModularData", "MergerData"});

    // df.Describe().Print();

    // Filter bad events
    auto dfFilter = df.Filter( // Check for heavier clusters than Li
                          [](ActRoot::MergerData& m)
                          {
                              if(m.fHeavy.fQave > 2000.)
                                  return false;
                              return true;
                          },
                          {"MergerData"})
                        .Filter( // Check if heavy reaches end of ACTAR (other way to mask heavier clusters,
                                 // but can delete good events with heavy particle bad reconstructed)
                            [](ActRoot::MergerData& m)
                            {
                                auto TL {m.fHeavy.fTL};
                                auto theta {m.fThetaHeavy * TMath::DegToRad()};
                                auto z_end {m.fRP.X() + TL * TMath::Cos(theta)};
                                if(z_end < 240)
                                    return false;
                                return true;
                            },
                            {"MergerData"})
                        .Filter( // Most of the times high charge deposit
                                 // is masked by rp or is in beam cluster
                            [&](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                            {
                                auto rp {m.fRP};
                                auto rp_y {rp.Y() / 2};
                                // Run for all clusters
                                int counter {};
                                for(auto& cluster : tpc.fClusters)
                                {
                                    auto voxels {cluster.GetRefToVoxels()};
                                    for(auto& v : voxels)
                                    {
                                        if(v.GetPosition().Y() > rp_y - 3 &&
                                           v.GetPosition().Y() < rp_y + 3) // aprox L1 exclusion zone
                                            if(v.GetCharge() > 3000.)
                                                counter++;
                                    }
                                }
                                if(counter > 4)
                                    return false;
                                return true;
                            },
                            {"MergerData", "TPCData"});
    // .Filter( // Check if light cluster go to the pad or cathode
    //     [&](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
    //     {
    //         auto& rp {m.fRP}; // This equals aprox 110 mm
    //         // Get last point of the track
    //         auto idx = m.fLightIdx;
    //         if(idx < 0)
    //             return false;
    //         // guard against out-of-range index (this was causing vector::at to throw)
    //         if(static_cast<std::size_t>(idx) >= tpc.fClusters.size())
    //             return false;
    //         auto& cluster = tpc.fClusters.at(idx);
    //         auto line = cluster.GetLine();
    //         cluster.SortAlongDir(line.GetDirection().Unit());
    //         auto& voxels = cluster.GetRefToVoxels();
    //         auto lastZ = voxels.back().GetPosition().Z() * driftFactor;
    //         // Z difference relative to the RP
    //         //  deltaZ > 0  -> track goes to larger Z (away from pad plane)
    //         //  deltaZ < 0  -> track goes to smaller Z (closer to pad plane)
    //         double deltaZ = lastZ - rp.Z();
    //         if(deltaZ < -100 || deltaZ > 140) // aprox pad plane and cathode position
    //         {
    //             return false;
    //             // std::cout << "Event rejected by light cluster Z position: deltaZ = " << deltaZ <<
    //             // " mm\n"; std::cout << "RP Z position: " << rp.Z() << " mm, last point Z position:
    //             // " << lastZ << " mm\n";
    //         }
    //         // std::cout << "Event accepted by light cluster Z position: deltaZ = " << deltaZ << "
    //         // mm\n"; std::cout << "RP Z position: " << rp.Z() << " mm, last point Z position: " <<
    //         // lastZ << " mm\n";
    //         return true;
    //     },
    //     {"MergerData", "TPCData"})

    auto dfZandQ = dfFilter.Define("zDrift",
                                   [&](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                                   {
                                       auto& rp {m.fRP};

                                       auto idx = m.fLightIdx;
                                       if(idx < 0)
                                           return -100.0f;

                                       if(static_cast<std::size_t>(idx) >= tpc.fClusters.size())
                                           return -100.0f;

                                       auto& cluster = tpc.fClusters.at(idx);
                                       auto& voxels = cluster.GetRefToVoxels();

                                       if(tpc.fRPs.empty())
                                           return -100.0f;
                                       auto rpVox = tpc.fRPs.front();

                                       float maxDist = -1.0;
                                       float zExtreme = 0.0;

                                       for(auto& v : voxels)
                                       {
                                           auto pos = v.GetPosition();

                                           float dx = pos.X() - rpVox.X();
                                           float dy = pos.Y() - rpVox.Y();
                                           float dz = pos.Z() - rpVox.Z();

                                           float dist = std::sqrt(dx * dx + dy * dy + dz * dz);

                                           if(dist > maxDist)
                                           {
                                               maxDist = dist;
                                               zExtreme = pos.Z();
                                           }
                                       }

                                       // diferencia en Z respecto al RP voxel
                                       float deltaZ = (zExtreme - rpVox.Z()) * driftFactor;

                                       float zDrift = 110.0 + deltaZ;

                                       return zDrift;
                                   },
                                   {"MergerData", "TPCData"});

    auto dfFilterZandQ =
        dfZandQ.Filter([](float zDrift, ActRoot::MergerData& m)
                       { return zDrift > -2.99441 && zDrift < 251.461 && m.fLight.fQave > 129.208; },
                       {"zDrift", "MergerData"}); // Keep only events with light cluster in the detector
    //////////////////////////////////////////////////////////
    // Get events in elastic area of plot theta vs Qtotal
    //////////////////////////////////////////////////////////

    // Get the cut for PID
    ActRoot::CutsManager<std::string> cuts;
    cuts.ReadCut("l1", TString::Format("./Cuts/%s_TLvsQ_%s.root", light.c_str(), beam.c_str()).Data());
    cuts.ReadCut("l1_theta", TString::Format("./Cuts/%s_ThetaVSq_%s.root", light.c_str(), beam.c_str()).Data());

    // Apply cut in theta - Qtotal
    auto dfDeuterium = dfFilter.Filter(
        [&](ActRoot::MergerData& m) { return cuts.IsInside("l1", m.fLight.fRawTL, m.fLight.fQtotal); }, {"MergerData"});
    auto dfElastic =
        dfDeuterium.Filter([&](ActRoot::MergerData& m)
                           { return cuts.IsInside("l1_theta", m.fThetaLight, m.fLight.fQtotal); }, {"MergerData"});

    auto dfPID = dfElastic.Define("IsDeuterium",
                                  [&](ActRoot::MergerData& m)
                                  {
                                      auto& hProfile {m.fQProf};
                                      // Fix bin errors to 1 to use chi2 test without statistical errors
                                      for(int i = 1; i <= hProfile.GetNbinsX(); ++i)
                                      {
                                          double content = hProfile.GetBinContent(i);
                                          double error = hProfile.GetBinError(i);
                                          hProfile.SetBinError(i, 1);
                                      }
                                      // double bestScore = 1e20;
                                      // double bestManualScore = 1e20;
                                      double bestFitScore = 1e20;
                                      // std::string bestParticle;
                                      // std::string bestManualParticle;
                                      std::string bestFitParticle;
                                      for(const auto& key : particles)
                                      {
                                          auto fSpline = FitSRIMtoChargeProfileFixedEnd(&hProfile, splineMap[key], key);

                                          if(!fSpline)
                                              continue;

                                          // ---------- ROOT fit chi2 ----------
                                          double chiFit = fSpline->GetChisquare() / fSpline->GetNDF();
                                          if(chiFit < bestFitScore)
                                          {
                                              bestFitScore = chiFit;
                                              bestFitParticle = key;
                                          }

                                          // // ---------- build model histogram ----------
                                          // TH1* model = BuildModelHistogramFromTF1(&hProfile, fSpline, "model_" +
                                          // key);
                                          //
                                          // // ---------- compare shape ----------
                                          // TH1* dataNorm = (TH1*)hProfile.Clone(("dataNorm_" + key).c_str());
                                          // TH1* modelNorm = (TH1*)model->Clone(("modelNorm_" + key).c_str());
                                          //
                                          // NormalizeHistogram(dataNorm);
                                          // NormalizeHistogram(modelNorm);
                                          //
                                          // // now is normalice inside function of chi2
                                          // double chiManually = Chi2Manually(&hProfile, model, fSpline);
                                          // double chiShape = Chi2Shape(&hProfile, model, fSpline);
                                          //
                                          // if(chiManually < bestManualScore)
                                          // {
                                          //     bestManualScore = chiManually;
                                          //     bestManualParticle = key;
                                          // }
                                          // if(chiShape < bestScore)
                                          // {
                                          //     bestScore = chiShape;
                                          //     bestParticle = key;
                                          // }
                                      }
                                      return bestFitParticle; // deuterium is the one we want to identify
                                  },
                                  {"MergerData"});

    /////////////////////////////////////////
    // Correct the range with the BP position
    /////////////////////////////////////////

    auto dfRange = dfPID.Define("Range",
                                [&](ActRoot::MergerData& m, std::string key)
                                {
                                    auto& hProfile {m.fQProf};
                                    // Fix bin errors to 1 to use chi2 test without statistical errors
                                    for(int i = 1; i <= hProfile.GetNbinsX(); ++i)
                                    {
                                        double content = hProfile.GetBinContent(i);
                                        double error = hProfile.GetBinError(i);
                                        hProfile.SetBinError(i, 1);
                                    }
                                    TF1* bestFitFunction =
                                        FitSRIMtoChargeProfileFixedEnd(&hProfile, splineMap[key], key);
                                    // Get range from spline max - fitted end shift
                                    double range = 0;
                                    if(bestFitFunction)
                                    {
                                        // Get the maximum of the TF1 (Bragg peak)
                                        // Then get the value where Eloss decreases to a 1/5 of the maximum after it
                                        double max = bestFitFunction->GetMaximum();
                                        double maxPos = bestFitFunction->GetMaximumX();
                                        double rangeValue = max / 5.0;

                                        for(double x = maxPos; x < bestFitFunction->GetXmax(); x += 0.05)
                                        {
                                            double y = bestFitFunction->Eval(x);
                                            if(y < rangeValue)
                                            {
                                                range = x;
                                                break;
                                            }
                                        }
                                    }
                                    return range;
                                },
                                {"MergerData", "IsDeuterium"});


    ////////////////////////////////////////////////////////////
    // Get the proportion of the best fit and the deuterium fit
    ////////////////////////////////////////////////////////////
    auto dfChiNormalized =
        dfPID.Define("chi2NormalizedToDeuterium",
                     [&](ActRoot::MergerData& m, std::string pid)
                     {
                         auto& hProfile {m.fQProf};
                         // Fix bin errors to 1 to use chi2 test without statistical errors
                         for(int i = 1; i <= hProfile.GetNbinsX(); ++i)
                         {
                             double content = hProfile.GetBinContent(i);
                             double error = hProfile.GetBinError(i);
                             hProfile.SetBinError(i, 1);
                         }
                         auto fDeuterium = FitSRIMtoChargeProfileFixedEnd(&hProfile, splineMap["light"], "light");
                         auto fParticle = FitSRIMtoChargeProfileFixedEnd(&hProfile, splineMap[pid], pid);
                         if(!fDeuterium || !fParticle)
                             return 1e12;
                         double chiDeuterium = fDeuterium->GetChisquare() / fDeuterium->GetNDF();
                         double chiParticle = fParticle->GetChisquare() / fParticle->GetNDF();
                         return chiParticle / chiDeuterium;
                     },
                     {"MergerData", "IsDeuterium"});


    // Count how many deuterium we have
    auto countDeuterium = dfPID.Filter([](const std::string& pid) { return pid == "lightD"; }, {"IsDeuterium"}).Count();
    auto countTritium = dfPID.Filter([](const std::string& pid) { return pid == "lightT"; }, {"IsDeuterium"}).Count();
    auto countProton = dfPID.Filter([](const std::string& pid) { return pid == "light"; }, {"IsDeuterium"}).Count();
    // Get the proportion of deuterium
    auto totalCount = dfPID.Count();
    std::cout << "Deuterium count: " << *countDeuterium << "\n";
    std::cout << "Tritium count: " << *countTritium << "\n";
    std::cout << "Proton count: " << *countProton << "\n";
    std::cout << "Total count: " << *totalCount << "\n";

    // Plot PID with and without cuts
    auto hTL_Qtot = dfFilter.Histo2D(
        {"hTL_Qtot", "Track Length vs Qtotal;Track Length [a.u.];Qtotal [a.u.]", 240, 0, 120, 2000, 0, 3e5},
        "MergerData.fLight.fRawTL", "MergerData.fLight.fQtotal");
    auto hTL_Qtot_cut =
        dfElastic.Histo2D({"hTL_Qtot_cut", "Track Length vs Qtotal with cut;Track Length [a.u.];Qtotal [a.u.]", 240, 0,
                           120, 2000, 0, 3e5},
                          "MergerData.fLight.fRawTL", "MergerData.fLight.fQtotal");
    auto hTL_Qtot_filterZandQ = dfFilterZandQ.Histo2D(
        {"hTL_Qtot_filterZandQ", "Track Length vs Qtotal with Z and Q cut;Track Length [a.u.];Qtotal [a.u.]", 240, 0,
         120, 2000, 0, 3e5},
        "MergerData.fLight.fRawTL", "MergerData.fLight.fQtotal");

    auto hTheta_Qtot =
        dfFilter.Histo2D({"hTheta_Qtot", "Theta vs Qtotal;Theta [deg];Qtotal [a.u.]", 180, 0, 180, 2000, 0, 3e5},
                         "MergerData.fThetaLight", "MergerData.fLight.fQtotal");
    auto hTheta_Qtot_cut = dfElastic.Histo2D(
        {"hTheta_Qtot_cut", "Theta vs Qtotal with cut;Theta [deg];Qtotal [a.u.]", 180, 0, 180, 2000, 0, 3e5},
        "MergerData.fThetaLight", "MergerData.fLight.fQtotal");
    auto hTheta_Qtot_filterZandQ =
        dfFilterZandQ.Histo2D({"hTheta_Qtot_filterZandQ", "Theta vs Qtotal with Z and Q cut;Theta [deg];Qtotal [a.u.]",
                               180, 0, 180, 2000, 0, 3e5},
                              "MergerData.fThetaLight", "MergerData.fLight.fQtotal");

    auto hPID_correctedRange =
        dfRange.Histo2D({"hPID_correctedRange", "Corrected Range vs Qtotal;Corrected Range [mm];Qtotal [a.u.]", 120, 0,
                         240, 2000, 0, 3e5},
                        "Range", "MergerData.fLight.fQtotal");
    auto hPID_no_correctedRange = dfRange.Histo2D(
        {"hPID_norm_correctedRange", "Corrected Range vs Qtotal normalized;Corrected Range [mm];Qtotal [a.u.]", 120, 0,
         240, 2000, 0, 3e5},
        "MergerData.fLight.fTL", "MergerData.fLight.fQtotal");

    auto hChi2NormalizedByDeuterium = dfChiNormalized.Histo1D(
        {"hChi2NormalizedByDeuterium", "Chi2 normalized to proton Chi2;Chi2 normalized to proton;Counts", 300, 0, 3},
        "chi2NormalizedToDeuterium");

    auto* c1 = new TCanvas("c1", "c1", 1200, 800);
    c1->Divide(2, 2);
    c1->cd(1);
    hTL_Qtot->DrawClone("COLZ");
    cuts.DrawCut("l1");
    c1->cd(2);
    hTL_Qtot_cut->DrawClone("COLZ");
    c1->cd(3);
    hTheta_Qtot->DrawClone("COLZ");
    cuts.DrawCut("l1_theta");
    c1->cd(4);
    hTheta_Qtot_cut->DrawClone("COLZ");

    auto cChi2 = new TCanvas("cChi2", "cChi2", 800, 600);
    hChi2NormalizedByDeuterium->DrawClone();

    auto* c2 = new TCanvas("c2", "c2", 800, 600);
    c2->Divide(2, 1);
    c2->cd(1);
    hPID_correctedRange->DrawClone("COLZ");
    c2->cd(2);
    hPID_no_correctedRange->DrawClone("COLZ");

    auto* c3 = new TCanvas("Filter Z and Q comparison", "Filter Z and Q comparison", 800, 600);
    c3->Divide(2, 2);
    c3->cd(1);
    hTheta_Qtot->DrawClone("COLZ");
    c3->cd(2);
    hTheta_Qtot_filterZandQ->DrawClone("COLZ");
    c3->cd(3);
    hTL_Qtot->DrawClone("COLZ");
    c3->cd(4);
    hTL_Qtot_filterZandQ->DrawClone("COLZ");

    // auto gr_lightD = new TGraph();
    // auto gr_lightT = new TGraph();
    // auto gr_light = new TGraph();
    // gr_lightD->SetName("gr_lightD");
    // gr_lightT->SetName("gr_lightT");
    // gr_light->SetName("gr_light");
    // gr_lightD->SetMarkerColor(kRed);
    // gr_lightT->SetMarkerColor(kBlue);
    // gr_light->SetMarkerColor(kGreen + 2);
    // gr_lightD->SetMarkerStyle(6);
    // gr_lightT->SetMarkerStyle(6);
    // gr_light->SetMarkerStyle(6);
    //
    // dfPID.Foreach(
    //     [&](ActRoot::MergerData& m, const std::string& pid)
    //     {
    //         if(pid == "lightD")
    //             gr_lightD->SetPoint(gr_lightD->GetN(), m.fThetaLight, m.fLight.fQtotal);
    //         else if(pid == "lightT")
    //             gr_lightT->SetPoint(gr_lightT->GetN(), m.fThetaLight, m.fLight.fQtotal);
    //         else if(pid == "light")
    //             gr_light->SetPoint(gr_light->GetN(), m.fThetaLight, m.fLight.fQtotal);
    //     },
    //     {"MergerData", "IsDeuterium"});
    //
    // auto gr_lightD_PID = new TGraph();
    // auto gr_lightT_PID = new TGraph();
    // auto gr_light_PID = new TGraph();
    // gr_lightD_PID->SetName("gr_lightD_PID");
    // gr_lightT_PID->SetName("gr_lightT_PID");
    // gr_light_PID->SetName("gr_light_PID");
    // gr_lightD_PID->SetMarkerColor(kRed);
    // gr_lightT_PID->SetMarkerColor(kBlue);
    // gr_light_PID->SetMarkerColor(kGreen + 2);
    // gr_lightD_PID->SetMarkerStyle(6);
    // gr_lightT_PID->SetMarkerStyle(6);
    // gr_light_PID->SetMarkerStyle(6);
    //
    // dfPID.Foreach(
    //     [&](ActRoot::MergerData& m, const std::string& pid)
    //     {
    //         if(pid == "lightD")
    //             gr_lightD_PID->SetPoint(gr_lightD_PID->GetN(), m.fLight.fRawTL, m.fLight.fQtotal);
    //         else if(pid == "lightT")
    //             gr_lightT_PID->SetPoint(gr_lightT_PID->GetN(), m.fLight.fRawTL, m.fLight.fQtotal);
    //         else if(pid == "light")
    //             gr_light_PID->SetPoint(gr_light_PID->GetN(), m.fLight.fRawTL, m.fLight.fQtotal);
    //     },
    //     {"MergerData", "IsDeuterium"});
    //
    // auto* c3 = new TCanvas("c3", "c3", 800, 600);
    // c3->Divide(2, 2);
    // c3->cd(1);
    // gr_lightD->SetTitle("All particles 3 runs candidates;Theta [deg];Qtotal [a.u.]");
    // gr_lightD->Draw("AP");
    // gr_lightT->SetTitle("Tritium candidates;Theta [deg];Qtotal [a.u.]");
    // gr_lightT->Draw("P SAME");
    // gr_light->SetTitle("Proton candidates;Theta [deg];Qtotal [a.u.]");
    // gr_light->Draw("P SAME");
    // TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    // legend->AddEntry(gr_lightD, "Deuterium", "P");
    // legend->AddEntry(gr_lightT, "Tritium", "P");
    // legend->AddEntry(gr_light, "Proton", "P");
    // legend->Draw();
    // c3->cd(2);
    // gr_lightD->SetTitle("Deuterium candidates;Theta [deg];Qtotal [a.u.]");
    // gr_lightD->Draw("AP");
    // c3->cd(3);
    // gr_lightT->SetTitle("Tritium candidates;Theta [deg];Qtotal [a.u.]");
    // gr_lightT->Draw("AP");
    // c3->cd(4);
    // gr_light->SetTitle("Proton candidates;Theta [deg];Qtotal [a.u.]");
    // gr_light->Draw("AP");
    //
    // auto* c4 = new TCanvas("c4", "c4", 800, 600);
    // c4->Divide(2, 2);
    // c4->cd(1);
    // gr_lightD_PID->SetTitle("Deuterium candidates;Track Length [a.u.];Qtotal [a.u.]");
    // gr_lightD_PID->Draw("AP");
    // gr_lightT_PID->SetTitle("Tritium candidates;Track Length [a.u.];Qtotal [a.u.]");
    // gr_lightT_PID->Draw("P SAME");
    // gr_light_PID->SetTitle("Proton candidates;Track Length [a.u.];Qtotal [a.u.]");
    // gr_light_PID->Draw("P SAME");
    // TLegend* legendPID = new TLegend(0.7, 0.7, 0.9, 0.9);
    // legendPID->AddEntry(gr_lightD_PID, "Deuterium", "P");
    // legendPID->AddEntry(gr_lightT_PID, "Tritium", "P");
    // legendPID->AddEntry(gr_light_PID, "Proton", "P");
    // legendPID->Draw();
    // c4->cd(2);
    // gr_lightD_PID->SetTitle("Deuterium candidates;Track Length [a.u.];Qtotal [a.u.]");
    // gr_lightD_PID->Draw("AP");
    // c4->cd(3);
    // gr_lightT_PID->SetTitle("Tritium candidates;Track Length [a.u   .];Qtotal [a.u.]");
    // gr_lightT_PID->Draw("AP");
    // c4->cd(4);
    // gr_light_PID->SetTitle("Proton candidates;Track Length [a.u.];Qtotal [a.u.]");
    // gr_light_PID->Draw("AP");

    // Fit the tritiums in a ForEach to plot the fit result over the hProfile
    // auto cProfile = new TCanvas("cProfile", "cProfile", 800, 600);
    // int counter = 0;
    // dfPID.Foreach(
    //     [&](ActRoot::MergerData& m, const std::string& pid)
    //     {
    //         if(pid == "light")
    //         {
    //             auto& hProfile {m.fQProf};
    //             auto fSplineT = FitSRIMtoChargeProfileFixedEnd(&hProfile, splineMap["lightT"], "lightT");
    //             auto fSplineD = FitSRIMtoChargeProfileFixedEnd(&hProfile, splineMap["lightD"], "lightD");
    //             auto fSplineP = FitSRIMtoChargeProfileFixedEnd(&hProfile, splineMap["light"], "light");
    //             if(fSplineT)
    //             {
    //                 TCanvas* c = new TCanvas(("c_" + pid + "_" + std::to_string(counter)).c_str(),
    //                                          ("Fit for " + pid).c_str(), 800, 600);
    //                 hProfile.DrawClone("hist");
    //                 fSplineT->SetLineColor(kGreen);
    //                 fSplineT->DrawClone("same P");
    //                 fSplineD->SetLineColor(kBlack);
    //                 fSplineD->DrawClone("same");
    //                 fSplineP->SetLineColor(kBlue);
    //                 fSplineP->DrawClone("same");
    //
    //                 // Create text for chi2 values, normalize to lower chi2 to 1 for better visualization
    //                 double chi2T = fSplineT->GetChisquare() / fSplineT->GetNDF();
    //                 double chi2D = fSplineD->GetChisquare() / fSplineD->GetNDF();
    //                 double chi2P = fSplineP->GetChisquare() / fSplineP->GetNDF();
    //                 double minChi2 = std::min({chi2T, chi2D, chi2P});
    //                 chi2T /= minChi2;
    //                 chi2D /= minChi2;
    //                 chi2P /= minChi2;
    //                 TLatex* text =
    //                     new TLatex(0.15, 0.85, Form("#chi^{2}/NDF: T=%.4f, D=%.4f, P=%.4f", chi2T, chi2D, chi2P));
    //                 text->SetNDC();
    //                 text->Draw();
    //             }
    //             counter++;
    //         }
    //     },
    //     {"MergerData", "IsDeuterium"});

    // Save events with good deuterium fit to ofstream
    // std::ofstream outFileD("./Outputs/good_deuterium_events.dat");
    // std::ofstream outFileT("./Outputs/good_tritium_events.dat");
    // std::ofstream outFileP("./Outputs/good_proton_events.dat");
    // dfPID.Foreach(
    //     [&](ActRoot::MergerData& m, const std::string& pid)
    //     {
    //         if(pid == "lightD")
    //         {
    //             m.Stream(outFileD);
    //         }
    //         else if(pid == "lightT")
    //         {
    //             m.Stream(outFileT);
    //         }
    //         else if(pid == "light")
    //         {
    //             m.Stream(outFileP);
    //         }
    //     },
    //     {"MergerData", "IsDeuterium"});
}