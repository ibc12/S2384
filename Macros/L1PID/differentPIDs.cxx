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

    // Define some variables to compare different PIDs
    auto def1 = def.Define("chargePerVoxelLight",
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
                           {"MergerData", "TPCData"});
    auto def2 = def1.Define("qRatio",
                            [&](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                            {
                                auto idx = m.fLightIdx;
                                if(idx < 0)
                                    return -1.0f;
                                auto& cluster = tpc.fClusters[idx];
                                auto& voxels = cluster.GetRefToVoxels();

                                // Dirección del track desde la línea ajustada
                                auto dir = cluster.GetLine().GetDirection().Unit();
                                auto rp = tpc.fRPs.front(); // punto de reacción

                                float qHead = 0, qTail = 0;

                                for(auto& v : voxels)
                                {
                                    auto pos = v.GetPosition();
                                    float proj = (pos - rp).Dot(dir); // distancia proyectada desde RP
                                    float halfTL = m.fLight.fRawTL / 2.0f;
                                    if(proj < halfTL)
                                    {
                                        qHead += v.GetCharge();
                                    }
                                    else
                                    {
                                        qTail += v.GetCharge();
                                    }
                                }
                                if(qHead < 1.0f)
                                    return -1.0f;
                                return qTail / qHead; // NORMALIZAR POR nVOXELS en cada parte?
                            },
                            {"MergerData", "TPCData"});
    auto df = def2.Define("TLsquare", [&](ActRoot::MergerData& m) { return m.fLight.fRawTL * m.fLight.fRawTL; },
                          {"MergerData"});

    // Get the cut for PID
    ActRoot::CutsManager<std::string> cuts;
    cuts.ReadCut("l1", TString::Format("./Cuts/%s_TLvsQ_%s.root", light.c_str(), beam.c_str()).Data());
    cuts.ReadCut("l1_theta", TString::Format("./Cuts/%s_ThetaVSq_%s.root", light.c_str(), beam.c_str()).Data());
    cuts.ReadCut("l1_chargePerVoxel",
                 TString::Format("./Cuts/%s_TLvsQvoxels_%s.root", light.c_str(), beam.c_str()).Data());

    // Apply cut in theta - Qtotal
    auto dfLight = df.Filter([&](ActRoot::MergerData& m)
                             { return cuts.IsInside("l1", m.fLight.fRawTL, m.fLight.fQtotal); }, {"MergerData"});
    auto dfAngle = df.Filter([&](ActRoot::MergerData& m)
                             { return cuts.IsInside("l1_theta", m.fThetaLight, m.fLight.fQtotal); }, {"MergerData"});
    auto dfVoxelLightCut = dfLight.Filter([&](ActRoot::MergerData& m, float qVoxel)
                                          { return cuts.IsInside("l1_chargePerVoxel", m.fThetaLight, qVoxel); },
                                          {"MergerData", "chargePerVoxelLight"});

    // PIDs and angle histograms
    auto h_TLvsQtotal = df.Histo2D({"h_TLvsQtotal", "L1 PID;Raw TL [a.u.];Qtotal [a.u.]", 240, 0, 120, 2000, 0, 3e5},
                                   "MergerData.fLight.fRawTL", "MergerData.fLight.fQtotal");
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

    // Try different variable plots
    auto h_TLsquareVsQtotal =
        df.Histo2D({"h_TLsquareVsQtotal", "L1 PID;Raw TL^{2} [a.u.];Qtotal [a.u.]", 240, 0, 120 * 120, 2000, 0, 3e5},
                   "TLsquare", "MergerData.fLight.fQtotal");
    auto h_qRatioVsQtotal =
        df.Histo2D({"h_qRatioVsQtotal", "L1 PID;Q_{tail}/Q_{head};Qtotal [a.u.]", 240, 0, 10, 2000, 0, 3e5}, "qRatio",
                   "MergerData.fLight.fQtotal");
    auto h_qRatioVsTheta =
        df.Histo2D({"h_qRatioVsTheta", "L1 PID;#theta_{Light} [#circ];Q_{tail}/Q_{head}", 240, 0, 180, 240, 0, 10},
                   "MergerData.fThetaLight", "qRatio");
    auto h_TLVsqRatio =
        df.Histo2D({"h_TLVsqRatio", "L1 PID;Raw TL;Qratio [a.u.]", 240, 0, 120, 240, 0, 10},
                   "MergerData.fLight.fRawTL", "qRatio");


    // Plot
    auto* c1 = new TCanvas("Canvas Data and 1st cut", "Canvas data and 1st cut", 1200, 800);
    c1->Divide(3, 2);
    c1->cd(1);
    h_TLvsQtotal->DrawClone("colz");
    cuts.DrawCut("l1");
    c1->cd(2);
    h_TLvsQvoxels->DrawClone("colz");
    c1->cd(3);
    h_thetaVsQ->DrawClone("colz");
    c1->cd(4);
    h_TLvsQtotal_cut->DrawClone("colz");
    c1->cd(5);
    h_TLvsQvoxels_cut->DrawClone("colz");
    cuts.DrawCut("l1_chargePerVoxel");
    c1->cd(6);
    h_thetaVsQtotal_cut->DrawClone("colz");

    auto* c2 = new TCanvas("Canvas angle cut", "Canvas angle cut", 1200, 400);
    c2->Divide(3, 1);
    c2->cd(1);
    h_TLvsQtotal_cut2->DrawClone("colz");
    c2->cd(2);
    h_TLvsQvoxels_cut2->DrawClone("colz");
    cuts.DrawCut("l1_chargePerVoxel");
    c2->cd(3);
    h_thetaVsQtotal_cut2->DrawClone("colz");

    auto* c3 = new TCanvas("Canvas 3rd cut", "Canvas 3rd cut", 1200, 400);
    c3->Divide(3, 1);
    c3->cd(1);
    h_TLvsQtotal_cut3->DrawClone("colz");
    c3->cd(2);
    h_TLvsQvoxels_cut3->DrawClone("colz");
    // cuts.DrawCut("l1_chargePerVoxel");
    c3->cd(3);
    h_thetaVsQtotal_cut3->DrawClone("colz");

    auto* c4 = new TCanvas("Canvas other variables", "Canvas other variables", 1200, 400);
    c4->Divide(2, 2);
    c4->cd(1);
    h_TLsquareVsQtotal->DrawClone("colz");
    c4->cd(2);
    h_qRatioVsQtotal->DrawClone("colz");
    c4->cd(3);
    h_qRatioVsTheta->DrawClone("colz");
    c4->cd(4);
    h_TLVsqRatio->DrawClone("colz");
}