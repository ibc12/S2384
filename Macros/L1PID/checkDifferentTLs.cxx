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


void ScalePoint(ROOT::Math::XYZPointF& point, float xy, float z, bool addOffset)
{
    if(addOffset) // when converting a bin point to physical units which wasnt already corrected
        point += ROOT::Math::XYZVectorF {0.5, 0.5, 0.5};
    point.SetX(point.X() * xy);
    point.SetY(point.Y() * xy);
    point.SetZ(point.Z() * z);
}

void checkDifferentTLs()
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

    auto df =
        dforigin.Filter([](ActRoot::ModularData& mod, ActRoot::MergerData& mer)
                        { return mod.Get("GATCONF") == 8 && (mer.fLightIdx != -1); }, {"ModularData", "MergerData"});

    auto dfoldTL = df.Define("oldTL",
                             [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                             {
                                 ROOT::Math::XYZPointF begin {tpc.fRPs.front()};
                                 auto lightCl {tpc.fClusters[m.fLightIdx]};
                                 auto line {lightCl.GetLine()};
                                 // Sort voxels in case they are not
                                 lightCl.SortAlongDir();
                                 auto end {lightCl.GetVoxels().back().GetPosition()};
                                 ScalePoint(end, 1, 1, true); // Put the position of the voxel on its center
                                 // Get projections onto fit
                                 auto projBegin {line.ProjectionPointOnLine(begin)};
                                 auto projEnd {line.ProjectionPointOnLine(end)};
                                 // Get TL as the distance between the projections
                                 return (projEnd - projBegin).R();
                             },
                             {"MergerData", "TPCData"})
                       .Define("nVoxelsLight",
                               [&](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                               {
                                   auto idx = m.fLightIdx;
                                   if(idx < 0)
                                       return -1;
                                   int nVoxels = tpc.fClusters[idx].GetSizeOfVoxels();
                                   return nVoxels;
                               },
                               {"MergerData", "TPCData"})
                       .Define("distanceRPtoFirstVoxel",
                               [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                               {
                                   ROOT::Math::XYZPointF rp {tpc.fRPs.front()};
                                   auto lightCl {tpc.fClusters[m.fLightIdx]};
                                   auto line {lightCl.GetLine()};
                                   // Sort voxels in case they are not
                                   lightCl.SortAlongDir();
                                   auto firstVoxel {lightCl.GetVoxels().front().GetPosition()};
                                   ScalePoint(firstVoxel, 1, 1, true); // Put the position of the voxel on its center
                                   // Get projections onto fit
                                   auto projBegin {line.ProjectionPointOnLine(rp)};
                                   auto projEnd {line.ProjectionPointOnLine(firstVoxel)};
                                   // Return distance between the projections
                                   return (projEnd - projBegin).R();
                               },
                               {"MergerData", "TPCData"});
    auto dfOldTL_phiPositive = dfoldTL.Filter([](ActRoot::MergerData& m) { return m.fPhiLight > 0; }, {"MergerData"});
    auto dfOldTL_phiNegative = dfoldTL.Filter([](ActRoot::MergerData& m) { return m.fPhiLight < 0; }, {"MergerData"});

    // DO some cuts
    ActRoot::CutsManager<std::string> cuts;
    cuts.ReadCut("d_oldTL_Q_restricted", "./Cuts/d_OldTLvsQ_restricted_11Li.root"); // Area of two bananas
    auto dfOldTL_cut = dfoldTL.Filter([&](ActRoot::MergerData& m, float oldTL)
                                      { return cuts.IsInside("d_oldTL_Q_restricted", oldTL, m.fLight.fQtotal); },
                                      {"MergerData", "oldTL"});
    auto dfOldTL_cut_phiPositive =
        dfOldTL_cut.Filter([](ActRoot::MergerData& m) { return m.fPhiLight > 0; }, {"MergerData"});
    auto dfOldTL_cut_phiNegative =
        dfOldTL_cut.Filter([](ActRoot::MergerData& m) { return m.fPhiLight < 0; }, {"MergerData"});

    // Do the plots to see the PID differences
    auto h_TLvsQtotal = df.Histo2D({"h_TLvsQtotal", "L1 PID;Raw TL [a.u.];Qtotal [a.u.]", 240, 0, 120, 4000, 0, 3e5},
                                   "MergerData.fLight.fRawTL", "MergerData.fLight.fQtotal");
    auto h_oldTLvsQtotal =
        dfoldTL.Histo2D({"h_oldTLvsQtotal", "L1 PID with old TL;Old TL [mm];Qtotal [a.u.]", 240, 0, 120, 4000, 0, 3e5},
                        "oldTL", "MergerData.fLight.fQtotal");
    // Plots of both TL to try to get deendances
    auto h_TLvsoldTL = dfoldTL.Histo2D({"h_TLvsoldTL", "L1 PID;Raw TL [a.u.];Old TL [mm]", 240, 0, 120, 240, 0, 120},
                                       "MergerData.fLight.fRawTL", "oldTL");
    auto h_TLvsoldTL_phiPositive = dfOldTL_phiPositive.Histo2D(
        {"h_TLvsoldTL_phiPositive", "L1 PID for phi positive;Raw TL [a.u.];Old TL [mm]", 240, 0, 120, 240, 0, 120},
        "MergerData.fLight.fRawTL", "oldTL");
    auto h_TLvsoldTL_phiNegative = dfOldTL_phiNegative.Histo2D(
        {"h_TLvsoldTL_phiNegative", "L1 PID for phi negative;Raw TL [a.u.];Old TL [mm]", 240, 0, 120, 240, 0, 120},
        "MergerData.fLight.fRawTL", "oldTL");
    // Plots of TLs individually to see diference
    auto h_oldTL = dfOldTL_cut.Histo1D({"h_oldTL", "L1 PID with old TL;Old TL [mm]", 240, 0, 120}, "oldTL");
    auto h_oldTL_phiPositive = dfOldTL_cut_phiPositive.Histo1D(
        {"h_oldTL_phiPositive", "L1 PID with old TL for phi positive;Old TL [mm];Counts", 240, 0, 120}, "oldTL");
    auto h_oldTL_phiNegative = dfOldTL_cut_phiNegative.Histo1D(
        {"h_oldTL_phiNegative", "L1 PID with old TL for phi negative;Old TL [mm];Counts", 240, 0, 120}, "oldTL");
    auto h_RawTL = dfOldTL_cut.Histo1D({"h_RawTL", "L1 PID with Raw TL;Raw TL [a.u.];Counts", 240, 0, 120},
                                       "MergerData.fLight.fRawTL");
    auto h_RawTL_phiPositive = dfOldTL_cut_phiPositive.Histo1D(
        {"h_RawTL_phiPositive", "L1 PID with Raw TL for phi positive;Raw TL [a.u.];Counts", 240, 0, 120},
        "MergerData.fLight.fRawTL");
    auto h_RawTL_phiNegative = dfOldTL_cut_phiNegative.Histo1D(
        {"h_RawTL_phiNegative", "L1 PID with Raw TL for phi negative;Raw TL [a.u.];Counts", 240, 0, 120},
        "MergerData.fLight.fRawTL");
    // Plots of number of voxels
    auto h_nVoxels = dfOldTL_cut.Histo1D(
        {"h_nVoxels", "Number of voxels in the light cluster;N voxels;Counts", 100, 0, 100}, "nVoxelsLight");
    auto h_nVoxels_phiPositive = dfOldTL_cut_phiPositive.Histo1D(
        {"h_nVoxels_phiPositive", "Number of voxels in the light cluster for phi positive;N voxels;Counts", 100, 0,
         100},
        "nVoxelsLight");
    auto h_nVoxels_phiNegative = dfOldTL_cut_phiNegative.Histo1D(
        {"h_nVoxels_phiNegative", "Number of voxels in the light cluster for phi negative;N voxels;Counts", 100, 0,
         100},
        "nVoxelsLight");
    // Plots of distance between both begins in the TLs definitions
    auto h_distanceRPtoFirstVoxel = dfOldTL_cut.Histo1D(
        {"h_distanceRPtoFirstVoxel", "Distance between RP and first voxel;Distance [mm];Counts", 100, 0, 20},
        "distanceRPtoFirstVoxel");
    auto h_distanceRPtoFirstVoxel_phiPositive = dfOldTL_cut_phiPositive.Histo1D(
        {"h_distanceRPtoFirstVoxel_phiPositive",
         "Distance between RP and first voxel for phi positive;Distance [mm];Counts", 100, 0, 20},
        "distanceRPtoFirstVoxel");
    auto h_distanceRPtoFirstVoxel_phiNegative = dfOldTL_cut_phiNegative.Histo1D(
        {"h_distanceRPtoFirstVoxel_phiNegative",
         "Distance between RP and first voxel for phi negative;Distance [mm];Counts", 100, 0, 20},
        "distanceRPtoFirstVoxel");


    auto* c = new TCanvas("cPIDcomparison", "Comparison of PID with different TL definitions", 1200, 400);
    c->Divide(2, 1);
    c->cd(1);
    h_TLvsQtotal->DrawClone("colz");
    c->cd(2);
    h_oldTLvsQtotal->DrawClone("colz");

    auto* c2 = new TCanvas("cTLcomparison", "Comparison of TL definitions", 1200, 400);
    c2->Divide(3, 1);
    c2->cd(1);
    h_TLvsoldTL->DrawClone("colz");
    c2->cd(2);
    h_TLvsoldTL_phiPositive->DrawClone("colz");
    c2->cd(3);
    h_TLvsoldTL_phiNegative->DrawClone("colz");

    auto* c3 = new TCanvas("cTLdistributions", "Comparison of TL distributions", 1200, 400);
    c3->Divide(3, 2);
    c3->cd(1);
    h_oldTL->DrawClone("colz");
    c3->cd(2);
    h_oldTL_phiPositive->DrawClone("colz");
    c3->cd(3);
    h_oldTL_phiNegative->DrawClone("colz");
    c3->cd(4);
    h_RawTL->DrawClone("colz");
    c3->cd(5);
    h_RawTL_phiPositive->DrawClone("colz");
    c3->cd(6);
    h_RawTL_phiNegative->DrawClone("colz");

    auto* c4 = new TCanvas("cNVoxelsAndDistance", "Number of voxels and distance between TL definitions", 1200, 400);
    c4->Divide(3, 2);
    c4->cd(1);
    h_nVoxels->DrawClone("colz");
    c4->cd(2);
    h_nVoxels_phiPositive->DrawClone("colz");
    c4->cd(3);
    h_nVoxels_phiNegative->DrawClone("colz");
    c4->cd(4);
    h_distanceRPtoFirstVoxel->DrawClone("colz");
    c4->cd(5);
    h_distanceRPtoFirstVoxel_phiPositive->DrawClone("colz");
    c4->cd(6);
    h_distanceRPtoFirstVoxel_phiNegative->DrawClone("colz");

    // Save some events
    // ActRoot::CutsManager<std::string> cuts;
    // cuts.ReadCut("RawTL_TLold", "./Cuts/RawTL_bigger_OldTL.root");
    // std::ofstream outFile("./Outputs/events_RawTLbiggerOldTL.dat");
    // dfOldTL_phiNegative.Foreach(
    //     [&](ActRoot::MergerData& m, float oldTL)
    //     {
    //         if(cuts.IsInside("RawTL_TLold", m.fLight.fRawTL, oldTL))
    //         {
    //             m.Stream(outFile);
    //         }
    //     },
    //     {"MergerData", "oldTL"});
}