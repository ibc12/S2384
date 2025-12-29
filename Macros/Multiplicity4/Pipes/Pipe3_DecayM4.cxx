#ifndef Pipe3_Deca7M4_cxx
#define Pipe3_Deca7M4_cxx

#include "ActCluster.h"
#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActParticle.h"
#include "ActSRIM.h"
#include "ActSilData.h"
#include "ActSilMatrix.h"
#include "ActSilSpecs.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TH2.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TString.h"

#include <map>
#include <string>

#include "../../../PostAnalysis/HistConfig.h"
#include "../../../PrettyStyle.C"
#include "../Utils.h"

struct DecayInfo
{
    int maxThetaIdx;
    int minThetaIdx;
    int maxQLengthIdx;
    int minQLengthIdx;

    double maxTheta;
    double minTheta;
    double maxQLength;
    double minQLength;

    double beamQLength;
    double lightQLength;
};

double ComputeQLength(ActRoot::Cluster& cl)
{
    cl.SortAlongDir(cl.GetLine().GetDirection());
    auto voxels = cl.GetRefToVoxels();

    if(voxels.size() < 2)
        return -1.0;

    float totalCharge = 0.f;
    for(const auto& v : voxels)
        totalCharge += v.GetCharge();

    ROOT::Math::XYZPointF firstPos = voxels.front().GetPosition();
    ROOT::Math::XYZPointF lastPos = voxels.back().GetPosition();
    Utils::ScalePoint(firstPos);
    Utils::ScalePoint(lastPos);

    double length = (lastPos - firstPos).R();
    if(length <= 0.0)
        return -1.0;

    return totalCharge / length;
}

void Pipe3_DecayM4(const std::string& beam, const std::string& target, const std::string& light)
{
    // Get file from pipe2
    TString infile = TString::Format("./Outputs/ExM4_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str());
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df("Final_Tree", infile.Data());

    // Create branch isnide Final_Tree for decay information


    // Analysis of the two decay particles
    auto dfDecay = df.Define("Decay",
                             [](int lightIdx, int beamIdx, ActRoot::TPCData& tpc)
                             {
                                 DecayInfo info;
                                 auto& clusters = tpc.fClusters;

                                 // Inicializar valores por defecto
                                 info.maxThetaIdx = -1;
                                 info.minThetaIdx = -1;
                                 info.maxQLengthIdx = -1;
                                 info.minQLengthIdx = -1;

                                 info.maxTheta = -1.0;
                                 info.minTheta = 1e9;
                                 info.maxQLength = -1.0;
                                 info.minQLength = 1e9;

                                 info.beamQLength = -1.0;
                                 info.lightQLength = -1.0;

                                 // Dirección del beam
                                 auto beamDir = clusters[beamIdx].GetLine().GetDirection().Unit();

                                 // ---- Beam Q/L ----
                                 info.beamQLength = ComputeQLength(clusters[beamIdx]);

                                 // ---- Light Q/L ----
                                 info.lightQLength = ComputeQLength(clusters[lightIdx]);

                                 // ---- Loop sobre otros clusters para decays ----
                                 for(int i = 0; i < clusters.size(); ++i)
                                 {
                                     if(i == lightIdx || i == beamIdx)
                                         continue;

                                     // ---- Angular info ----
                                     auto currentDir = clusters[i].GetLine().GetDirection().Unit();
                                     double angle = TMath::ACos(beamDir.Dot(currentDir)) * TMath::RadToDeg();

                                     if(angle > info.maxTheta)
                                     {
                                         info.maxTheta = angle;
                                         info.maxThetaIdx = i;
                                     }
                                     if(angle < info.minTheta)
                                     {
                                         info.minTheta = angle;
                                         info.minThetaIdx = i;
                                     }

                                     // ---- Q/L info ----
                                     double ql = ComputeQLength(clusters[i]);
                                     if(ql < 0.0)
                                         continue;

                                     if(ql > info.maxQLength)
                                     {
                                         info.maxQLength = ql;
                                         info.maxQLengthIdx = i;
                                     }
                                     if(ql < info.minQLength)
                                     {
                                         info.minQLength = ql;
                                         info.minQLengthIdx = i;
                                     }
                                 }

                                 return info;
                             },
                             {"LightIdx", "BeamIdx", "TPCData"});

    // Try to get silicon information for the decay particles
    

    auto dfPlot = dfDecay.Define("beamQLength", "Decay.beamQLength")
                      .Define("lightQLength", "Decay.lightQLength")
                      .Define("maxQLength", "Decay.maxQLength")
                      .Define("minQLength", "Decay.minQLength");

    // Build histos
    auto hChargeBeam = dfPlot.Histo1D({"hChargeBeam", "Charge/Distance;Counts", 150, 0, 2000}, "beamQLength");
    auto hChargeLight = dfPlot.Histo1D({"hChargeLight", "Charge/Distance;Counts", 150, 0, 2000}, "lightQLength");
    auto hChargeMax = dfPlot.Histo1D({"hChargeMax", "Charge/Distance;Counts", 150, 0, 2000}, "maxQLength");
    auto hChargeMin = dfPlot.Histo1D({"hChargeMin", "Charge/Distance;Counts", 150, 0, 2000}, "minQLength");

    // Plot and legend
    auto c0 = new TCanvas("cDecayM4_0", "Decay Q/L comparison", 800, 600);
    c0->cd();

    hChargeBeam->SetLineColor(kBlue + 1);
    hChargeLight->SetLineColor(kMagenta + 1);
    hChargeMax->SetLineColor(kRed + 1);
    hChargeMin->SetLineColor(kGreen + 2);

    auto hBeamDrawn = (TH1D*)hChargeBeam->DrawClone();
    auto hLightDrawn = (TH1D*)hChargeLight->DrawClone("same");
    auto hMaxDrawn = (TH1D*)hChargeMax->DrawClone("same");
    auto hMinDrawn = (TH1D*)hChargeMin->DrawClone("same");

    auto leg = new TLegend(0.60, 0.65, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    leg->AddEntry(hBeamDrawn, "Beam Q/L", "l");
    leg->AddEntry(hLightDrawn, "Light Q/L", "l");
    leg->AddEntry(hMaxDrawn, "Decay max Q/L", "l");
    leg->AddEntry(hMinDrawn, "Decay min Q/L", "l");
    leg->Draw();


    // Save dataframe in a .root file
    TString outfile = TString::Format("./Outputs/DecayM4_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str());
    dfDecay.Snapshot("Final_Tree", outfile);
    std::cout << "Saving Final_Tree in " << outfile << "from Pipe2_ExM4 " << '\n';
}
#endif