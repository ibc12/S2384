#include "ActCluster.h"
#include "ActContinuity.h"
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


struct UnrebinedVoxel
{
    UnrebinedVoxel(int nPadsX, int nPadsY, int nPadsZ) : nPadsX(nPadsX), nPadsY(nPadsY), nPadsZ(nPadsZ) {}

    std::map<unsigned int, std::vector<ActRoot::Voxel>>
        rebinnedIndexToVoxels {}; // Coorelation of Global Index in rebinned space to voxels in original space
    int nPadsX {};
    int nPadsY {};
    int nPadsZ {};

    unsigned int BuildGlobalIndex(const int& x, const int& y, const int& z)
    {
        return x + y * nPadsX + z * nPadsX * nPadsY;
    }
};


std::pair<std::vector<ActRoot::Voxel>, UnrebinedVoxel> RebinTracks(const std::vector<ActRoot::Voxel>& voxels,
                                                                   ActRoot::TPCParameters* tpcPars, int rebinX,
                                                                   int rebinY = -1, int rebinZ = -1)
{
    // Rebin all voxels by the rebin factor and store the rebined voxels (no repeat positions, if so sum the charges)
    std::vector<ActRoot::Voxel> rebinnedVoxels;
    int nPadsX = tpcPars->GetNPADSX() / rebinX;
    if(rebinY == -1)
        rebinY = rebinX;
    if(rebinZ == -1)
        rebinZ = rebinX;
    int nPadsY = tpcPars->GetNPADSY() / rebinY;
    int nPadsZ = tpcPars->GetNPADSZ() / rebinZ;
    UnrebinedVoxel unrebinedData {nPadsX, nPadsY, nPadsZ};
    // std::cout << "pad numbers: " << unrebinedData.nPadsX << " " << unrebinedData.nPadsY << " " << unrebinedData.nPadsZ
    //           << '\n';
    // position
    std::map<unsigned int, unsigned int> rebinedIndexAndPosition; // global rebined index -> position in rebinnedVoxels
                                                                  // vector
    for(auto& voxel : voxels)
    {
        auto pos = voxel.GetPosition();
        int x = std::floor(pos.X() / rebinX);
        int y = std::floor(pos.Y() / rebinY);
        int z = std::floor(pos.Z() / rebinZ);

        // Get the global index of the rebinned voxel
        auto globalIndex = unrebinedData.BuildGlobalIndex(x, y, z);
        unrebinedData.rebinnedIndexToVoxels[globalIndex].push_back(voxel);

        if(!rebinedIndexAndPosition.count(globalIndex))
        {
            ActRoot::Voxel rebinnedVoxel {{static_cast<float>(x), static_cast<float>(y), static_cast<float>(z)},
                                          voxel.GetCharge(),
                                          voxel.GetIsSaturated()};
            rebinnedVoxels.push_back(rebinnedVoxel);
            rebinedIndexAndPosition[globalIndex] = rebinnedVoxels.size() - 1;
        }
        else
        {
            auto idx {rebinedIndexAndPosition[globalIndex]};
            rebinnedVoxels[idx].SetCharge(rebinnedVoxels[idx].GetCharge() + voxel.GetCharge());
        }
        // std::cout << "rebinnedVoxels size: " << rebinnedVoxels.size() << '\n';
    }
    return {rebinnedVoxels, unrebinedData};
}

std::vector<ActRoot::Voxel> UndoRebinning(std::vector<ActRoot::Voxel>& rebinnedVoxels, UnrebinedVoxel& rebinnedData)
{
    // Compute the global index for the rebined voxel, and return the voxels in the map with that index
    std::vector<ActRoot::Voxel> unrebinedVoxels;
    for(const auto& voxel : rebinnedVoxels)
    {
        auto pos = voxel.GetPosition();
        int x = static_cast<int>(pos.X());
        int y = static_cast<int>(pos.Y());
        int z = static_cast<int>(pos.Z());

        unsigned int globalIndex = rebinnedData.BuildGlobalIndex(x, y, z);
        if(rebinnedData.rebinnedIndexToVoxels.count(globalIndex))
        {
            unrebinedVoxels.insert(unrebinedVoxels.end(), rebinnedData.rebinnedIndexToVoxels.at(globalIndex).begin(),
                                   rebinnedData.rebinnedIndexToVoxels.at(globalIndex).end());
        }
    }
    return unrebinedVoxels;
}

void fixBrokenTracksZ()
{
    // Get one event that has this kind of problem r66e36114
    std::string dataconf {"../../configs/data.conf"};
    ActRoot::DataManager dataman {dataconf, ActRoot::ModeType::EMerge};
    dataman.SetRuns(66, 66);
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    auto chain3 {dataman.GetChain(ActRoot::ModeType::EFilter)};
    auto chain4 {dataman.GetChain(ActRoot::ModeType::EReadTPC)};
    chain->AddFriend(chain2.get());
    chain->AddFriend(chain3.get(), "TPCData");
    chain->AddFriend(chain4.get(), "GETTree");

    ROOT::RDataFrame df {*chain};
    // df.Describe().Print();

    // Get the parameters for TPCDetector
    ActRoot::TPCParameters* tpcPars = new ActRoot::TPCParameters("Actar");
    std::cout << "TPC Parameters: " << '\n';
    tpcPars->SetREBINZ(4);
    tpcPars->Print();

    int rebinFactor = 4;
    ActRoot::TPCParameters tpcParsRebined = {
        (tpcPars->GetNPADSX() / rebinFactor),
        (tpcPars->GetNPADSY() / rebinFactor),
        (tpcPars->GetNPADSZ() / rebinFactor),
    };
    tpcParsRebined.Print();

    ActAlgorithm::Continuity continuity {&tpcParsRebined, 2};

    auto dfEvent = df.Filter([](const ActRoot::MergerData& m) { return m.fEntry == 36114; }, {"MergerData"});

    // Print angles of light to see if that is the event
    dfEvent.Foreach(
        [](const ActRoot::MergerData& m)
        {
            std::cout << "ThetaLight: " << m.fThetaLight << '\n';
            std::cout << "PhiLight: " << m.fPhiLight << '\n';
        },
        {"MergerData"});

    // Get the TPCData of that event and print the raw voxels
    dfEvent.Foreach([](const ActRoot::TPCData& tpc) { std::cout << "Raw voxels: " << tpc.fRaw.size() << '\n'; },
                    {"GETTree.TPCData"});

    auto* hEventXY = new TH2D("hEventXY", "Event display;X;Y", tpcPars->GetNPADSX(), 0, tpcPars->GetNPADSX(),
                              tpcPars->GetNPADSY(), 0, tpcPars->GetNPADSY());
    auto* hNoiseXY = new TH2D("hNoiseXY", "Noise display;X;Y", tpcPars->GetNPADSX(), 0, tpcPars->GetNPADSX(),
                              tpcPars->GetNPADSY(), 0, tpcPars->GetNPADSY());
    auto* hAllVoxelsXY = new TH2D("hAllVoxelsXY", "All voxels display;X;Y", tpcPars->GetNPADSX(), 0,
                                  tpcPars->GetNPADSX(), tpcPars->GetNPADSY(), 0, tpcPars->GetNPADSY());

    auto* hEventXZ = new TH2D("hEventXZ", "Event display;X;Z", tpcPars->GetNPADSX(), 0, tpcPars->GetNPADSX(),
                              tpcPars->GetNPADSZ(), 0, tpcPars->GetNPADSZ());
    auto* hNoiseXZ = new TH2D("hNoiseXZ", "Noise display;X;Z", tpcPars->GetNPADSX(), 0, tpcPars->GetNPADSX(),
                              tpcPars->GetNPADSZ(), 0, tpcPars->GetNPADSZ());
    auto* hAllVoxelsXZ = new TH2D("hAllVoxelsXZ", "All voxels display;X;Z", tpcPars->GetNPADSX(), 0,
                                  tpcPars->GetNPADSX(), tpcPars->GetNPADSZ(), 0, tpcPars->GetNPADSZ());

    auto* hEventYZ = new TH2D("hEventYZ", "Event display;Y;Z", tpcPars->GetNPADSY(), 0, tpcPars->GetNPADSY(),
                              tpcPars->GetNPADSZ(), 0, tpcPars->GetNPADSZ());
    auto* hNoiseYZ = new TH2D("hNoiseYZ", "Noise display;Y;Z", tpcPars->GetNPADSY(), 0, tpcPars->GetNPADSY(),
                              tpcPars->GetNPADSZ(), 0, tpcPars->GetNPADSZ());
    auto* hAllVoxelsYZ = new TH2D("hAllVoxelsYZ", "All voxels display;Y;Z", tpcPars->GetNPADSY(), 0,
                                  tpcPars->GetNPADSY(), tpcPars->GetNPADSZ(), 0, tpcPars->GetNPADSZ());

    auto* hEventAfterProcessingXY =
        new TH2D("hEventAfterProcessingXY", "Event display after processing;X;Y", tpcPars->GetNPADSX(), 0,
                 tpcPars->GetNPADSX(), tpcPars->GetNPADSY(), 0, tpcPars->GetNPADSY());
    auto* hEventAfterProcessingXZ =
        new TH2D("hEventAfterProcessingXZ", "Event display after processing;X;Z", tpcPars->GetNPADSX(), 0,
                 tpcPars->GetNPADSX(), tpcPars->GetNPADSZ(), 0, tpcPars->GetNPADSZ());
    auto* hEventAfterProcessingYZ =
        new TH2D("hEventAfterProcessingYZ", "Event display after processing;Y;Z", tpcPars->GetNPADSY(), 0,
                 tpcPars->GetNPADSY(), tpcPars->GetNPADSZ(), 0, tpcPars->GetNPADSZ());

    std::vector<ActRoot::Voxel> voxelsBefore;
    std::vector<ActRoot::Voxel> voxelsAfter;

    // Now, rebin the TPC data in Z, apply continuity, and then undo the rebinning to get the corrected raw voxels.
    dfEvent.Foreach(
        [&](ActRoot::TPCData& tpcFilter, ActRoot::TPCData& tpc, ActRoot::MergerData& m)
        {
            // Get the raw voxels and the light cluster
            int lightIdx = m.fLightIdx;
            if(lightIdx == -1)
            {
                std::cout << "No light cluster found, skipping event" << '\n';
                return;
            }
            auto lightCluster = tpcFilter.fClusters[lightIdx];
            std::cout << "Light cluster before processing: " << lightCluster.GetVoxels().size() << " voxels" << '\n';
            lightCluster.GetLine().Print();
            auto lightVoxels = lightCluster.GetVoxels();
            auto rawVoxels = tpc.fRaw;
            std::cout << "Light voxels: " << lightVoxels.size() << '\n';
            std::cout << "Raw voxels: " << rawVoxels.size() << '\n';
            // Get a std::vector of clusters with light and raw voxels
            std::vector<ActRoot::Voxel> voxels;
            voxels.insert(voxels.end(), lightVoxels.begin(), lightVoxels.end());
            voxels.insert(voxels.end(), rawVoxels.begin(), rawVoxels.end());

            // Rebin in Z
            auto result = RebinTracks(voxels, tpcPars, rebinFactor, rebinFactor, rebinFactor);

            // Apply continuity
            auto clusterRebinedAfterContinuity = continuity.Run(result.first, true);
            // Get the vector of voxels of the cluster, if there is more than one cluster, get the voxels of the biggest
            // one
            std::vector<ActRoot::Voxel> VoxelsToUnrebin;
            if(!clusterRebinedAfterContinuity.first.empty())
            {
                auto maxCluster = std::max_element(clusterRebinedAfterContinuity.first.begin(),
                                                   clusterRebinedAfterContinuity.first.end(),
                                                   [](const ActRoot::Cluster& a, const ActRoot::Cluster& b)
                                                   { return a.GetVoxels().size() < b.GetVoxels().size(); });
                VoxelsToUnrebin = maxCluster->GetRefToVoxels();
            }
            else
            {
                std::cout << "No clusters found after continuity, skipping event" << '\n';
                return;
            }

            // Undo rebinning
            auto unrebinedClusterAfterContinuity = UndoRebinning(VoxelsToUnrebin, result.second);

            // Get clusters to do the plots
            ActRoot::Cluster clusterAfterProcessing {};
            clusterAfterProcessing.SetVoxels(unrebinedClusterAfterContinuity);
            clusterAfterProcessing.SortAlongDir(lightCluster.GetLine().GetDirection());
            clusterAfterProcessing.ReFit();
            std::cout << "Cluster after processing: " << clusterAfterProcessing.GetVoxels().size() << " voxels" << '\n';
            std::cout << "Clusters procesed: " << clusterRebinedAfterContinuity.first.size() << " voxels" << '\n';
            std::cout << "Noise after processing: " << clusterRebinedAfterContinuity.second.size() << " voxels" << '\n';
            clusterAfterProcessing.GetLine();

            for(auto& voxel : rawVoxels)
            {
                auto pos = voxel.GetPosition();
                hNoiseXY->Fill(pos.X(), pos.Y());
                hAllVoxelsXY->Fill(pos.X(), pos.Y());
                hNoiseXZ->Fill(pos.X(), pos.Z());
                hAllVoxelsXZ->Fill(pos.X(), pos.Z());
                hNoiseYZ->Fill(pos.Y(), pos.Z());
                hAllVoxelsYZ->Fill(pos.Y(), pos.Z());
            }
            for(auto& voxel : lightVoxels)
            {
                auto pos = voxel.GetPosition();
                hEventXY->Fill(pos.X(), pos.Y());
                hAllVoxelsXY->Fill(pos.X(), pos.Y());
                hEventXZ->Fill(pos.X(), pos.Z());
                hAllVoxelsXZ->Fill(pos.X(), pos.Z());
                hEventYZ->Fill(pos.Y(), pos.Z());
                hAllVoxelsYZ->Fill(pos.Y(), pos.Z());
            }
            for(auto& voxel : unrebinedClusterAfterContinuity)
            {
                auto pos = voxel.GetPosition();
                hEventAfterProcessingXY->Fill(pos.X(), pos.Y());
                hEventAfterProcessingXZ->Fill(pos.X(), pos.Z());
                hEventAfterProcessingYZ->Fill(pos.Y(), pos.Z());
            }

            // Set voxles before and after processing for later use
            voxelsBefore = lightVoxels;
            voxelsAfter = unrebinedClusterAfterContinuity;
        },
        {"TPCData", "GETTree.TPCData", "MergerData"});

    // Create clusters and get lines
    ActRoot::Cluster clusterBeforeProcessing {};
    clusterBeforeProcessing.SetVoxels(voxelsBefore);
    clusterBeforeProcessing.ReFit();
    ActRoot::Cluster clusterAfterProcessing {};
    clusterAfterProcessing.SetVoxels(voxelsAfter);
    clusterAfterProcessing.ReFit();

    // Get the lines
    auto lineBefore = clusterBeforeProcessing.GetLine();
    auto lineAfter = clusterAfterProcessing.GetLine();

    auto* c1 = new TCanvas("c1", "Event display", 1200, 400);
    c1->Divide(3, 1);

    // ===================== XY =====================
    c1->cd(1);

    // Cluster (azul)
    hEventXY->SetFillColorAlpha(kBlue, 0.35);
    hEventXY->SetLineColor(kBlue);

    // Noise (rojo)
    hNoiseXY->SetFillColorAlpha(kRed, 0.35);
    hNoiseXY->SetLineColor(kRed);

    // Dibujar
    hEventXY->Draw("BOX");
    hNoiseXY->Draw("BOX SAME");
    lineBefore.GetPolyLine("xy", 0, tpcPars->GetNPADSX(), tpcPars->GetNPADSY(), tpcPars->GetNPADSZ(), 1)
        ->DrawClone("same");

    // Leyenda
    {
        auto* leg = new TLegend(0.7, 0.75, 0.9, 0.9);
        leg->AddEntry(hEventXY, "Cluster", "f");
        leg->AddEntry(hNoiseXY, "Noise", "f");
        leg->Draw();
    }


    // ===================== XZ =====================
    c1->cd(2);

    hEventXZ->SetFillColorAlpha(kBlue, 0.35);
    hEventXZ->SetLineColor(kBlue);
    hNoiseXZ->SetFillColorAlpha(kRed, 0.35);
    hNoiseXZ->SetLineColor(kRed);

    hEventXZ->Draw("BOX");
    hNoiseXZ->Draw("BOX SAME");
    lineBefore.GetPolyLine("xz", 0, tpcPars->GetNPADSX(), tpcPars->GetNPADSY(), tpcPars->GetNPADSZ(), 1)
        ->DrawClone("same");


    {
        auto* leg = new TLegend(0.7, 0.75, 0.9, 0.9);
        leg->AddEntry(hEventXZ, "Cluster", "f");
        leg->AddEntry(hNoiseXZ, "Noise", "f");
        leg->Draw();
    }


    // ===================== YZ =====================
    c1->cd(3);

    hEventYZ->SetFillColorAlpha(kBlue, 0.35);
    hEventYZ->SetLineColor(kBlue);

    hNoiseYZ->SetFillColorAlpha(kRed, 0.35);
    hNoiseYZ->SetLineColor(kRed);

    hEventYZ->Draw("BOX");
    hNoiseYZ->Draw("BOX SAME");
    lineBefore.GetPolyLine("yz", 0, tpcPars->GetNPADSX(), tpcPars->GetNPADSY(), tpcPars->GetNPADSZ(), 1)
        ->DrawClone("same");


    {
        auto* leg = new TLegend(0.7, 0.75, 0.9, 0.9);
        leg->AddEntry(hEventYZ, "Cluster", "f");
        leg->AddEntry(hNoiseYZ, "Noise", "f");
        leg->Draw();
    }

    auto* c2 = new TCanvas("c2", "Event display after processing", 1200, 400);
    c2->Divide(3, 1);
    c2->cd(1);
    hEventAfterProcessingXY->Draw("BOX");
    auto polyLineBeforeXY =
        lineBefore.GetPolyLine("xy", 0, tpcPars->GetNPADSX(), tpcPars->GetNPADSY(), tpcPars->GetNPADSZ(), 1);
    polyLineBeforeXY->SetLineColor(kBlue);
    polyLineBeforeXY->DrawClone("same");
    lineAfter.GetPolyLine("xy", 0, tpcPars->GetNPADSX(), tpcPars->GetNPADSY(), tpcPars->GetNPADSZ(), 1)
        ->DrawClone("same");
    c2->cd(2);
    hEventAfterProcessingXZ->Draw("BOX");
    auto polyLineBeforeXZ =
        lineBefore.GetPolyLine("xz", 0, tpcPars->GetNPADSX(), tpcPars->GetNPADSY(), tpcPars->GetNPADSZ(), 1);
    polyLineBeforeXZ->SetLineColor(kBlue);
    polyLineBeforeXZ->DrawClone("same");
    lineAfter.GetPolyLine("xz", 0, tpcPars->GetNPADSX(), tpcPars->GetNPADSY(), tpcPars->GetNPADSZ(), 1)
        ->DrawClone("same");
    c2->cd(3);
    hEventAfterProcessingYZ->Draw("BOX");
    auto polyLineBeforeYZ =
        lineBefore.GetPolyLine("yz", 0, tpcPars->GetNPADSX(), tpcPars->GetNPADSY(), tpcPars->GetNPADSZ(), 1);
    polyLineBeforeYZ->SetLineColor(kBlue);
    polyLineBeforeYZ->DrawClone("same");
    lineAfter.GetPolyLine("yz", 0, tpcPars->GetNPADSX(), tpcPars->GetNPADSY(), tpcPars->GetNPADSZ(), 1)
        ->DrawClone("same");
}