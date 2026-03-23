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

unsigned int BuildGlobalIndex(const int& x, const int& y, const int& z, const ActRoot::TPCParameters& fPars)
{
    return x + y * fPars.GetNPADSX() + z * fPars.GetNPADSX() * fPars.GetNPADSY();
}

std::pair<std::vector<ActRoot::Cluster>, std::vector<ActRoot::Voxel>>
RebinTracks(std::vector<ActRoot::Voxel>& voxels, ActRoot::TPCParameters& fPars, int rebinFactor = 2)
{
    // Rebin in Z by a factor of rebinFactor
    for(auto& voxel : voxels)
    {
        auto pos = voxel.GetPosition();
        pos.SetX(std::floor(pos.X() / rebinFactor));
        pos.SetY(std::floor(pos.Y() / rebinFactor));
        pos.SetZ(std::floor(pos.Z() / rebinFactor));
        voxel.SetPosition(pos);
    }

    ActRoot::TPCParameters rebinnedPars {fPars.GetNPADSX() / rebinFactor, fPars.GetNPADSY() / rebinFactor,
                                         fPars.GetNPADSZ() / rebinFactor};

    // Apply continuity to the rebinned voxels.
    auto* continuity = new ActAlgorithm::Continuity(&rebinnedPars, 4);
    auto result = continuity->Run(voxels, true);
    std::cout << "Continuity result for noise: " << result.second.size() << " voxels" << '\n';
    std::cout << "Continuity result for clusters: " << result.first.size() << " clusters" << '\n';

    return result;
}

void UndoRebinning() {}

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
    df.Describe().Print();

    // Get the parameters for TPCDetector
    ActRoot::TPCParameters* tpcPars = new ActRoot::TPCParameters("Actar");
    std::cout << "TPC Parameters: " << '\n';
    tpcPars->Print();

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
                              tpcPars->GetNPADSZ() / 4, 0, tpcPars->GetNPADSZ() / 4);
    auto* hNoiseXZ = new TH2D("hNoiseXZ", "Noise display;X;Z", tpcPars->GetNPADSX(), 0, tpcPars->GetNPADSX(),
                              tpcPars->GetNPADSZ() / 4, 0, tpcPars->GetNPADSZ() / 4);
    auto* hAllVoxelsXZ = new TH2D("hAllVoxelsXZ", "All voxels display;X;Z", tpcPars->GetNPADSX(), 0,
                                  tpcPars->GetNPADSX(), tpcPars->GetNPADSZ() / 4, 0, tpcPars->GetNPADSZ() / 4);

    auto* hEventYZ = new TH2D("hEventYZ", "Event display;Y;Z", tpcPars->GetNPADSY(), 0, tpcPars->GetNPADSY(),
                              tpcPars->GetNPADSZ() / 4, 0, tpcPars->GetNPADSZ() / 4);
    auto* hNoiseYZ = new TH2D("hNoiseYZ", "Noise display;Y;Z", tpcPars->GetNPADSY(), 0, tpcPars->GetNPADSY(),
                              tpcPars->GetNPADSZ() / 4, 0, tpcPars->GetNPADSZ() / 4);
    auto* hAllVoxelsYZ = new TH2D("hAllVoxelsYZ", "All voxels display;Y;Z", tpcPars->GetNPADSY(), 0,
                                  tpcPars->GetNPADSY(), tpcPars->GetNPADSZ() / 4, 0, tpcPars->GetNPADSZ() / 4);

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
            auto lightVoxels = lightCluster.GetVoxels();
            auto rawVoxels = tpc.fRaw;
            std::cout << "Light voxels: " << lightVoxels.size() << '\n';
            std::cout << "Raw voxels: " << rawVoxels.size() << '\n';
            // Get a std::vector of clusters with light and raw voxels
            std::vector<ActRoot::Voxel> voxels;
            voxels.insert(voxels.end(), lightVoxels.begin(), lightVoxels.end());
            voxels.insert(voxels.end(), rawVoxels.begin(), rawVoxels.end());

            // Rebin in Z
            auto result = RebinTracks(voxels, *tpcPars);
            for(auto& voxel : result.second)
            {
                auto pos = voxel.GetPosition();
                hNoiseXY->Fill(pos.X(), pos.Y());
                hAllVoxelsXY->Fill(pos.X(), pos.Y());
                hNoiseXZ->Fill(pos.X(), pos.Z());
                hAllVoxelsXZ->Fill(pos.X(), pos.Z());
                hNoiseYZ->Fill(pos.Y(), pos.Z());
                hAllVoxelsYZ->Fill(pos.Y(), pos.Z());
            }
            for(auto& cluster : result.first)
            {
                for(auto& voxel : cluster.GetVoxels())
                {
                    auto pos = voxel.GetPosition();
                    hEventXY->Fill(pos.X(), pos.Y());
                    hAllVoxelsXY->Fill(pos.X(), pos.Y());
                    hEventXZ->Fill(pos.X(), pos.Z());
                    hAllVoxelsXZ->Fill(pos.X(), pos.Z());
                    hEventYZ->Fill(pos.Y(), pos.Z());
                    hAllVoxelsYZ->Fill(pos.Y(), pos.Z());
                }
            }

            // Apply continuity

            // Undo rebinning
            UndoRebinning();
        },
        {"TPCData", "GETTree.TPCData", "MergerData"});

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

    {
        auto* leg = new TLegend(0.7, 0.75, 0.9, 0.9);
        leg->AddEntry(hEventYZ, "Cluster", "f");
        leg->AddEntry(hNoiseYZ, "Noise", "f");
        leg->Draw();
    }
}