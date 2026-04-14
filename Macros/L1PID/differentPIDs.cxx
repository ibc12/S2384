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


void differentPIDs()
{
    PrettyStyle(true, true);

    // Get all L1 experimental data for 7Li; cut in L1 and filter bad events, select elastic, check good fits to
    // deuterium

    // Get all data for 7Li (there is no triton there, so easier to see deuterium)
    std::string beam {"11Li"};
    std::string target {"d"};
    std::string light {"p"};

    auto inFile {TString::Format("../../PostAnalysis/Outputs/tree_preprocess_F_%s.root", beam.c_str())};

    // RDataFrame
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame dforigin {"PreProcessed_Tree", inFile};

    // Filter GATCONF == L1
    auto def =
        dforigin.Filter([](ActRoot::ModularData& mod, ActRoot::MergerData& mer)
                        { return mod.Get("GATCONF") == 8 && (mer.fLightIdx != -1); }, {"ModularData", "MergerData"});
    auto defSilicon =
        dforigin.Filter([](ActRoot::ModularData& mod) { return mod.Get("GATCONF") == 1; }, {"ModularData"})
            .Filter([](ActRoot::MergerData& mer) { return mer.fLightIdx != -1; }, {"MergerData"});

    // Define some variables to compare different PIDs
    auto dfDefinitions =
        def.Define("chargePerVoxelLight",
                   [&](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                   {
                       auto idx = m.fLightIdx;
                       if(idx < 0)
                           return -1.0f;
                       auto& voxels = tpc.fClusters[idx].GetRefToVoxels();
                       float totalCharge = 0.0f;
                       for(auto& v : voxels)
                           totalCharge += v.GetCharge();
                       return totalCharge / voxels.size();
                   },
                   {"MergerData", "TPCData"})
            .Define("qRatio",
                    [&](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                    {
                        auto idx = m.fLightIdx;
                        if(idx < 0)
                            return -1.0f;
                        auto& cluster = tpc.fClusters[idx];
                        auto& voxels = cluster.GetRefToVoxels();

                        // Track direction from the fitted line
                        auto dir = cluster.GetLine().GetDirection().Unit();
                        auto rp = tpc.fRPs.front(); // reaction point

                        float qHead = 0, qTail = 0;
                        int nHead = 0, nTail = 0;

                        for(auto& v : voxels)
                        {
                            auto pos = v.GetPosition();
                            float proj = (pos - rp).Dot(dir); // projected distance from RP
                            float halfTL = m.fLight.fRawTL / 2.0f;
                            if(proj < halfTL)
                            {
                                qHead += v.GetCharge();
                                nHead++;
                            }
                            else
                            {
                                qTail += v.GetCharge();
                                nTail++;
                            }
                        }
                        if(qHead < 1.0f)
                            return -1.0f;
                        return (qTail) / (qHead); // TODO: normalize by number of voxels in each part?
                    },
                    {"MergerData", "TPCData"})
            .Define("qRatioNorm",
                    [&](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                    {
                        auto idx = m.fLightIdx;
                        if(idx < 0)
                            return -1.0f;
                        auto& cluster = tpc.fClusters[idx];
                        auto& voxels = cluster.GetRefToVoxels();

                        // Track direction from the fitted line
                        auto dir = cluster.GetLine().GetDirection().Unit();
                        auto rp = tpc.fRPs.front(); // reaction point

                        float qHead = 0, qTail = 0;
                        int nHead = 0, nTail = 0;

                        for(auto& v : voxels)
                        {
                            auto pos = v.GetPosition();
                            float proj = (pos - rp).Dot(dir); // projected distance from RP
                            float halfTL = m.fLight.fRawTL / 2.0f;
                            if(proj < halfTL)
                            {
                                qHead += v.GetCharge();
                                nHead++;
                            }
                            else
                            {
                                qTail += v.GetCharge();
                                nTail++;
                            }
                        }
                        if(qHead < 1.0f)
                            return -1.0f;
                        return (qTail / nTail) / (qHead / nHead); // TODO: normalize by number of voxels in each part?
                    },
                    {"MergerData", "TPCData"})
            .Define("lastThirdQ",
                    [&](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                    {
                        auto idx = m.fLightIdx;
                        if(idx < 0)
                            return -1.0f;
                        auto& cluster = tpc.fClusters[idx];
                        auto& voxels = cluster.GetRefToVoxels();

                        // Track direction from the fitted line
                        auto dir = cluster.GetLine().GetDirection().Unit();
                        auto rp = tpc.fRPs.front(); // reaction point

                        float qLastThird = 0;
                        int nLastThird = 0;

                        for(auto& v : voxels)
                        {
                            auto pos = v.GetPosition();
                            float proj = (pos - rp).Dot(dir); // projected distance from RP
                            float twoThirdTL = 2.0f * m.fLight.fRawTL / 3.0f;
                            if(proj >= twoThirdTL)
                            {
                                qLastThird += v.GetCharge();
                                nLastThird++;
                            }
                        }
                        if(nLastThird == 0)
                            return -1.0f;
                        return qLastThird / nLastThird; // TODO: normalize by number of voxels in the last third?
                    },
                    {"MergerData", "TPCData"})
            .Define("TLsquare", [&](ActRoot::MergerData& m) { return m.fLight.fRawTL * m.fLight.fRawTL; },
                    {"MergerData"})
            .Define("nVoxelsLight",
                    [&](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                    {
                        auto idx = m.fLightIdx;
                        if(idx < 0)
                            return -1;
                        int nVoxels = tpc.fClusters[idx].GetSizeOfVoxels();
                        return nVoxels;
                    },
                    {"MergerData", "TPCData"});

    // Get the cut for PID
    ActRoot::CutsManager<std::string> cuts;
    cuts.ReadCut("l1_alfas", TString::Format("./Cuts/alfa_TL2vsQ_%s.root", beam.c_str()).Data());
    cuts.ReadCut("l1", TString::Format("./Cuts/%s_TLvsQ_%s.root", light.c_str(), beam.c_str()).Data());
    cuts.ReadCut("l1_theta", TString::Format("./Cuts/%s_ThetaVSq_%s.root", light.c_str(), beam.c_str()).Data());
    cuts.ReadCut("l1_chargePerVoxel",
                 TString::Format("./Cuts/%s_TLvsQvoxels_%s.root", light.c_str(), beam.c_str()).Data());
    // Cut for getting events
    cuts.ReadCut("heavyHighQave", "./Cuts/events_qAveHeavyHigh.root");
    cuts.ReadCut("phiFlat_lowTL", "./Cuts/events_phiFlat_LowTL.root");

    auto df = dfDefinitions.Filter([&](ActRoot::MergerData& m, float tlSquare)
                                   { return !cuts.IsInside("l1_alfas", tlSquare, m.fLight.fQtotal); },
                                   {"MergerData", "TLsquare"});

    // Apply cut in theta - Qtotal
    auto dfLight = df.Filter([&](ActRoot::MergerData& m)
                             { return cuts.IsInside("l1", m.fLight.fRawTL, m.fLight.fQtotal); }, {"MergerData"});
    auto dfAngle = df.Filter([&](ActRoot::MergerData& m)
                             { return cuts.IsInside("l1_theta", m.fThetaLight, m.fLight.fQtotal); }, {"MergerData"});
    auto dfVoxelLightCut = dfLight.Filter([&](ActRoot::MergerData& m, float qVoxel)
                                          { return cuts.IsInside("l1_chargePerVoxel", m.fThetaLight, qVoxel); },
                                          {"MergerData", "chargePerVoxelLight"});

    // Filter in phi to see the differences
    auto dfPhiPositive = df.Filter([&](ActRoot::MergerData& m) { return m.fPhiLight > 0; }, {"MergerData"});
    auto dfPhiNegative = df.Filter([&](ActRoot::MergerData& m) { return m.fPhiLight < 0; }, {"MergerData"});
    auto dfPhiFlat =
        df.Filter([&](ActRoot::MergerData& m)
                  { return (m.fPhiLight > 80 && m.fPhiLight < 100) || (m.fPhiLight > -100 && m.fPhiLight < -80); },
                  {"MergerData"});
    auto dfPhiVerticalPad =
        df.Filter([&](ActRoot::MergerData& m) { return (m.fPhiLight > 160 && m.fPhiLight < 180); }, {"MergerData"});
    auto dfPhiVerticalCathode =
        df.Filter([&](ActRoot::MergerData& m) { return (m.fPhiLight > -20 && m.fPhiLight < 0); }, {"MergerData"});


    // PIDs and angle histograms
    auto h_TLvsQtotal = df.Histo2D({"h_TLvsQtotal", "L1 PID;Raw TL [a.u.];Qtotal [a.u.]", 240, 0, 120, 2000, 0, 3e5},
                                   "MergerData.fLight.fRawTL", "MergerData.fLight.fQtotal");
    auto h_TLvsQave = df.Histo2D({"h_TLvsQave", "L1 PID;Raw TL [a.u.];Qave [a.u.]", 240, 0, 120, 2000, 0, 3500},
                                 "MergerData.fLight.fRawTL", "MergerData.fLight.fQave");
    auto h_TLvsQvoxels =
        df.Histo2D({"h_TLvsQvoxels", "L1 PID;Raw TL [a.u.];Charge per voxel [a.u.]", 240, 0, 120, 800, 0, 5000},
                   "MergerData.fLight.fRawTL", "chargePerVoxelLight");
    auto h_thetaVsQ =
        df.Histo2D({"h_thetaVsQ", "L1 #theta vs Q;#theta_{Light} [#circ];Qtotal [a.u.]", 240, 0, 180, 2000, 0, 3e5},
                   "MergerData.fThetaLight", "MergerData.fLight.fQtotal");

    // Apply first cut, in TLvsQtotal
    auto h_TLvsQtotal_cut = dfLight.Histo2D(
        {"h_TLvsQtotal_cut", "L1 PID cut in TLvsQ;Raw TL [a.u.];Qtotal [a.u.]", 240, 0, 120, 2000, 0, 3e5},
        "MergerData.fLight.fRawTL", "MergerData.fLight.fQtotal");
    auto h_TLvsQave_cut =
        dfLight.Histo2D({"h_TLvsQave_cut", "L1 PID cut in TLvsQ;Raw TL [a.u.];Qave [a.u.]", 240, 0, 120, 2000, 0, 3500},
                        "MergerData.fLight.fRawTL", "MergerData.fLight.fQave");
    auto h_TLvsQvoxels_cut =
        dfLight.Histo2D({"h_TLvsQvoxels_cut", "L1 PID cut in TLvsQvoxels;Raw TL [a.u.];Charge per voxel [a.u.]", 240, 0,
                         120, 800, 0, 5000},
                        "MergerData.fLight.fRawTL", "chargePerVoxelLight");
    auto h_thetaVsQtotal_cut =
        dfLight.Histo2D({"h_thetaVsQtotal_cut", "L1 #theta vs Q cut in TLvsQ;#theta_{Light} [#circ];Qtotal [a.u.]", 240,
                         0, 180, 2000, 0, 3e5},
                        "MergerData.fThetaLight", "MergerData.fLight.fQtotal");
    // Apply second cut, in TLvsQvoxels
    auto h_TLvsQtotal_cut2 = dfVoxelLightCut.Histo2D(
        {"h_TLvsQtotal_cut2", "L1 PID cut in TLvsQvoxels;Raw TL [a.u.];Qtotal [a.u.]", 240, 0, 120, 2000, 0, 3e5},
        "MergerData.fLight.fRawTL", "MergerData.fLight.fQtotal");
    auto h_TLvsQave_cut2 = dfVoxelLightCut.Histo2D(
        {"h_TLvsQave_cut2", "L1 PID cut in TLvsQvoxels;Raw TL [a.u.];Qave [a.u.]", 240, 0, 120, 2000, 0, 3500},
        "MergerData.fLight.fRawTL", "MergerData.fLight.fQave");
    auto h_TLvsQvoxels_cut2 = dfVoxelLightCut.Histo2D(
        {"h_TLvsQvoxels_cut2", "L1 PID cut in TLvsQvoxels;Raw TL [a.u.];Charge per voxel [a.u.]", 240, 0, 120, 200, 0,
         5000},
        "MergerData.fLight.fRawTL", "chargePerVoxelLight");
    auto h_thetaVsQtotal_cut2 = dfVoxelLightCut.Histo2D(
        {"h_thetaVsQtotal_cut2", "L1 #theta vs Q cut in TLvsQvoxels;#theta_{Light} [#circ];Qtotal [a.u.]", 240, 0, 180,
         2000, 0, 3e5},
        "MergerData.fThetaLight", "MergerData.fLight.fQtotal");
    // Apply cut in theta - Qtotal
    auto h_TLvsQtotal_cut3 = dfAngle.Histo2D(
        {"h_TLvsQtotal_cut3", "L1 PID cut in TLvsQvoxels;Raw TL [a.u.];Qtotal [a.u.]", 240, 0, 120, 2000, 0, 3e5},
        "MergerData.fLight.fRawTL", "MergerData.fLight.fQtotal");
    auto h_TLvsQvoxels_cut3 =
        dfAngle.Histo2D({"h_TLvsQvoxels_cut3", "L1 PID cut in TLvsQvoxels;Raw TL [a.u.];Charge per voxel [a.u.]", 240,
                         0, 120, 200, 0, 5000},
                        "MergerData.fLight.fRawTL", "chargePerVoxelLight");
    auto h_thetaVsQtotal_cut3 = dfAngle.Histo2D(
        {"h_thetaVsQtotal_cut3", "L1 #theta vs Q cut in TLvsQvoxels;#theta_{Light} [#circ];Qtotal [a.u.]", 240, 0, 180,
         2000, 0, 3e5},
        "MergerData.fThetaLight", "MergerData.fLight.fQtotal");
    auto h_TLvsQave_cut3 = dfAngle.Histo2D(
        {"h_TLvsQave_cut3", "L1 PID cut in TLvsQvoxels;Raw TL [a.u.];Qave [a.u.]", 240, 0, 120, 2000, 0, 3500},
        "MergerData.fLight.fRawTL", "MergerData.fLight.fQave");

    // Try different variable plots
    auto h_TLsquareVsQtotal =
        df.Histo2D({"h_TLsquareVsQtotal", "L1 PID;Raw TL^{2} [a.u.];Qtotal [a.u.]", 240, 0, 150 * 150, 4000, 0, 3e5},
                   "TLsquare", "MergerData.fLight.fQtotal");
    auto h_QLastThirdVsQtotal = df.Histo2D(
        {"h_QLastThirdVsQtotal", "L1 PID;Q in last third of track [a.u.];Qtotal [a.u.]", 2000, 0, 2e5, 2000, 0, 3e5},
        "lastThirdQ", "MergerData.fLight.fQtotal");
    auto h_qRatioVsQtotal =
        df.Histo2D({"h_qRatioVsQtotal", "L1 PID;Q_{tail}/Q_{head};Qtotal [a.u.]", 240, 0, 10, 2000, 0, 3e5}, "qRatio",
                   "MergerData.fLight.fQtotal");
    auto h_qRatioVsTheta =
        df.Histo2D({"h_qRatioVsTheta", "L1 PID;#theta_{Light} [#circ];Q_{tail}/Q_{head}", 240, 0, 180, 240, 0, 10},
                   "MergerData.fThetaLight", "qRatio");
    auto h_TLVsqRatio = df.Histo2D({"h_TLVsqRatio", "L1 PID;Raw TL;Qratio [a.u.]", 240, 0, 120, 240, 0, 10},
                                   "MergerData.fLight.fRawTL", "qRatio");
    auto h_qRatioNormVsQtotal =
        df.Histo2D({"h_qRatioNormVsQtotal", "L1 PID;Q_{tail}/Q_{head} Normaliced to voxels;Qtotal [a.u.]", 240, 0, 10,
                    2000, 0, 3e5},
                   "qRatioNorm", "MergerData.fLight.fQtotal");
    auto h_qRatioNormVsTheta =
        df.Histo2D({"h_qRatioNormVsTheta", "L1 PID;#theta_{Light} [#circ];Q_{tail}/Q_{head} Normaliced to voxels", 240,
                    0, 180, 240, 0, 10},
                   "MergerData.fThetaLight", "qRatioNorm");
    auto h_TLVsqRatioNorm = df.Histo2D(
        {"h_TLVsqRatioNorm", "L1 PID;Raw TL;Q_{tail}/Q_{head} Normaliced to voxels", 240, 0, 120, 240, 0, 10},
        "MergerData.fLight.fRawTL", "qRatioNorm");

    // Heavy particle correlations
    auto h_ThetaLight_ThetaHeavy = df.Histo2D(
        {"h_ThetaLight_ThetaHeavy", "Theta light vs Theta heavy for L1;#theta_{Light} [#circ];#theta_{Heavy} [#circ]",
         300, 0, 160, 90, 0, 30},
        "MergerData.fThetaLight", "MergerData.fThetaHeavy");
    auto h_QaveLight_QaveHeavy =
        df.Histo2D({"h_QaveLight_QaveHeavy", "Qave light vs Qave heavy for L1;Qave_{Light} [a.u.];Qave_{Heavy} [a.u.]",
                    1200, 0, 2500, 2000, 0, 3500},
                   "MergerData.fLight.fQave", "MergerData.fHeavy.fQave");
    auto h_QtotalLight_ThetaHeavy = df.Histo2D(
        {"h_QtotalLight_ThetaHeavy", "Qtotal light vs Theta heavy for L1;Qtotal_{Light} [a.u.];#theta_{Heavy} [#circ]",
         2000, 0, 3e5, 50, 0, 30},
        "MergerData.fLight.fQtotal", "MergerData.fThetaHeavy");

    // Phi different cases
    auto h_TLvsQtotal_PhiPositive = dfPhiPositive.Histo2D(
        {"h_TLvsQtotal_PhiPositive", "L1 PID phi > 0;Raw TL [a.u.];Qtotal [a.u.]", 240, 0, 120, 4000, 0, 3e5},
        "MergerData.fLight.fRawTL", "MergerData.fLight.fQtotal");
    auto h_TLvsQtotal_PhiNegative = dfPhiNegative.Histo2D(
        {"h_TLvsQtotal_PhiNegative", "L1 PID phi < 0;Raw TL [a.u.];Qtotal [a.u.]", 240, 0, 120, 4000, 0, 3e5},
        "MergerData.fLight.fRawTL", "MergerData.fLight.fQtotal");
    auto h_phi_silicon =
        defSilicon.Histo1D({"h_phi_silicon", "Phi distribution for silicon left;Phi [deg];Counts", 360, -180, 180},
                           "MergerData.fPhiLight");
    auto h_phi_L1 =
        df.Histo1D({"h_phi_L1", "Phi distribution for L1;Phi [deg];Counts", 360, -180, 180}, "MergerData.fPhiLight");
    auto h_TLvsQtotal_PhiFlat = dfPhiFlat.Histo2D(
        {"h_TLvsQtotal_PhiFlat", "L1 PID phi flat;Raw TL [a.u.];Qtotal [a.u.]", 240, 0, 120, 4000, 0, 3e5},
        "MergerData.fLight.fRawTL", "MergerData.fLight.fQtotal");
    auto h_TLvsQtotal_PhiVertical =
        dfPhiVerticalPad.Histo2D({"h_TLvsQtotal_PhiVertical pad", "L1 PID phi vertical pad;Raw TL [a.u.];Qtotal [a.u.]",
                                  240, 0, 120, 4000, 0, 3e5},
                                 "MergerData.fLight.fRawTL", "MergerData.fLight.fQtotal");
    auto h_TLvsQtotal_PhiVerticalCathode = dfPhiVerticalCathode.Histo2D(
        {"h_TLvsQtotal_PhiVertical cathode", "L1 PID phi vertical cathode;Raw TL [a.u.];Qtotal [a.u.]", 240, 0, 120,
         4000, 0, 3e5},
        "MergerData.fLight.fRawTL", "MergerData.fLight.fQtotal");
    auto h_nVoxelsPhi = dfDefinitions.Histo2D(
        {"h_nVoxelsPhi", "Number of voxels vs phi;Phi [deg];N voxels", 360, -180, 180, 100, 0, 100},
        "MergerData.fPhiLight", "nVoxelsLight");


    // Plot
    auto* c1 = new TCanvas("Canvas Data and 1st cut", "Canvas data and 1st cut", 1200, 1200);
    c1->Divide(4, 2);
    c1->cd(1);
    h_TLvsQtotal->DrawClone("colz");
    cuts.DrawCut("l1");
    c1->cd(2);
    h_TLvsQvoxels->DrawClone("colz");
    c1->cd(3);
    h_thetaVsQ->DrawClone("colz");
    c1->cd(4);
    h_TLvsQave->DrawClone("colz");
    c1->cd(5);
    h_TLvsQtotal_cut->DrawClone("colz");
    c1->cd(6);
    h_TLvsQvoxels_cut->DrawClone("colz");
    cuts.DrawCut("l1_chargePerVoxel");
    c1->cd(7);
    h_thetaVsQtotal_cut->DrawClone("colz");
    c1->cd(8);
    h_TLvsQave_cut->DrawClone("colz");

    auto* c2 = new TCanvas("Canvas angle cut", "Canvas angle cut", 1200, 1200);
    c2->Divide(4, 1);
    c2->cd(1);
    h_TLvsQtotal_cut2->DrawClone("colz");
    c2->cd(2);
    h_TLvsQvoxels_cut2->DrawClone("colz");
    cuts.DrawCut("l1_chargePerVoxel");
    c2->cd(3);
    h_thetaVsQtotal_cut2->DrawClone("colz");
    c2->cd(4);
    h_TLvsQave_cut2->DrawClone("colz");

    auto* c3 = new TCanvas("Canvas 3rd cut", "Canvas 3rd cut", 1200, 1200);
    c3->Divide(4, 1);
    c3->cd(1);
    h_TLvsQtotal_cut3->DrawClone("colz");
    c3->cd(2);
    h_TLvsQvoxels_cut3->DrawClone("colz");
    // cuts.DrawCut("l1_chargePerVoxel");
    c3->cd(3);
    h_thetaVsQtotal_cut3->DrawClone("colz");
    c3->cd(4);
    h_TLvsQave_cut3->DrawClone("colz");

    auto* c4 = new TCanvas("Canvas other variables", "Canvas other variables", 1200, 1200);
    c4->Divide(3, 3);
    c4->cd(1);
    h_TLsquareVsQtotal->DrawClone("colz");
    c4->cd(2);
    h_QLastThirdVsQtotal->DrawClone("colz");
    c4->cd(4);
    h_qRatioVsQtotal->DrawClone("colz");
    c4->cd(5);
    h_qRatioVsTheta->DrawClone("colz");
    c4->cd(6);
    h_TLVsqRatio->DrawClone("colz");
    c4->cd(7);
    h_qRatioNormVsQtotal->DrawClone("colz");
    c4->cd(8);
    h_qRatioNormVsTheta->DrawClone("colz");
    c4->cd(9);
    h_TLVsqRatioNorm->DrawClone("colz");

    auto* c5 = new TCanvas("Canvas heavy-light correlations", "Canvas heavy-light correlations", 1200, 1200);
    c5->Divide(2, 2);
    c5->cd(1);
    h_ThetaLight_ThetaHeavy->DrawClone("colz");
    c5->cd(2);
    h_QaveLight_QaveHeavy->DrawClone("colz");
    c5->cd(3);
    h_QtotalLight_ThetaHeavy->DrawClone("colz");

    auto* c6 = new TCanvas("Canvas phi correlation", "Canvas phi correlation", 1200, 600);
    c6->Divide(3, 2);
    c6->cd(1);
    h_TLvsQtotal_PhiPositive->DrawClone("colz");
    c6->cd(2);
    h_TLvsQtotal_PhiNegative->DrawClone("colz");
    c6->cd(3);
    h_TLvsQtotal->DrawClone("colz");
    // c6->cd(4);
    // h_phi_silicon->DrawClone();
    // c6->cd(5);
    // h_phi_L1->DrawClone();
    c6->cd(4);
    h_TLvsQtotal_PhiFlat->DrawClone("colz");
    c6->cd(5);
    h_TLvsQtotal_PhiVertical->DrawClone("colz");
    // c6->cd(6);
    // h_TLvsQtotal_PhiVerticalCathode->DrawClone("colz");
    c6->cd(6);
    h_nVoxelsPhi->DrawClone("colz");

    // Save some events
    cuts.ReadCut("deuterium", "./Cuts/d_TLvsQ_11Li.root");
    // std::ofstream outFile("./Outputs/eventsElastic_phiNear90positive.dat");
    // df.Foreach(
    //     [&](ActRoot::MergerData& m)
    //     {
    //         if(cuts.IsInside("deuterium", m.fLight.fRawTL, m.fLight.fQtotal))
    //         {
    //             if((m.fPhiLight > 80 && m.fPhiLight < 100) && (m.fThetaLight > 85 && m.fThetaLight < 90))
    //                 m.Stream(outFile);
    //         }
    //     },
    //     {"MergerData"});
}