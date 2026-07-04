#include "ActContinuity.h"
#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActTPCData.h"
#include "ActTPCParameters.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"

#include <fstream>
#include <iostream>

void BrokenTracksZL1()
{
    // ROOT::EnableImplicitMT();

    // ROOT::RDataFrame df {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_7Li_d_d_filtered.root"};
    ROOT::RDataFrame df {"PreProcessed_Tree", "../../PostAnalysis/Outputs/tree_preprocess_F_11Li.root"};

    // ActRoot::DataManager dataman {"../../configs/data_7Li.conf", ActRoot::ModeType::EMerge};
    // auto chain {dataman.GetChain()};
    // Add friends if necessary
    // auto friend1 {dataman.GetChain(ActRoot::ModeType::EFilter)};
    // chain->AddFriend(friend1.get());
    // auto friend2 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    // chain->AddFriend(friend2.get());
    //
    // ROOT::RDataFrame df {*chain};
    // auto def {
    //     df.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"})}; // only L1
    auto def {df.Filter([](ActRoot::ModularData& m) { return m.Get("GATCONF") == 8; }, {"ModularData"})}; // only L1

    // Print the number of entries in the dataframe
    auto nEntries = static_cast<int>(def.Count().GetValue());
    std::cout << "Number of entries in the dataframe: " << nEntries << std::endl;

    // Get tpc pars
    ActRoot::TPCParameters tpc {"Actar"};
    std::cout << "TPC parameters: " << std::endl;
    std::cout << " XLength = " << tpc.GetNPADSX() << std::endl;
    std::cout << " YLength = " << tpc.GetNPADSY() << std::endl;
    std::cout << " ZLength = " << tpc.GetNPADSZ() << std::endl;
    std::cout << " XPadSize = " << tpc.GetPadSide() << std::endl;

    // Build continuity object
    auto cont = std::make_shared<ActAlgorithm::Continuity>(&tpc, 5);
    // Go through the light particle and filter uncontinuous tracks
    auto defClusters = def.Define("nClustersInLight",
                                  [&](ActRoot::TPCData& tpc, ActRoot::MergerData& m)
                                  {
                                      auto lightIdx = m.fLightIdx;
                                      auto lightCl = tpc.fClusters[lightIdx];
                                      auto voxels = lightCl.GetVoxels();
                                      // Run continuity algorithm
                                      auto results = cont->Run(voxels);
                                      int nClusters = static_cast<int>(results.first.size());
                                      return nClusters;
                                  },
                                  {"TPCData", "MergerData"});
    auto defClustersFiltered = defClusters.Filter(
        [&](int nClustersInLight)
        {
            if(nClustersInLight > 1)
                return true; // uncontinuous track
            else
                return false; // continuous track
        },
        {"nClustersInLight"});

    auto defClusters0Filtered = defClusters.Filter(
        [&](int nClustersInLight)
        {
            if(nClustersInLight == 0)
                return true; // uncontinuous track
            else
                return false; // continuous track
        },
        {"nClustersInLight"});

    // Plot the phi and theta distributions of the uncontinuous tracks
    auto hPhi = defClustersFiltered.Histo1D(
        {"hPhi", "#phi distribution of uncontinuous tracks;#phi [rad];Counts", 100, -180, 180}, "fPhiLight");
    auto hTheta = defClustersFiltered.Histo1D(
        {"hTheta", "#theta distribution of uncontinuous tracks;#theta [rad];Counts", 100, 0, 180}, "fThetaLight");
    auto hPhiAll = defClusters.Histo1D({"hPhiAll", "#phi distribution of all tracks;#phi [rad];Counts", 100, -180, 180},
                                       "fPhiLight");
    auto hThetaAll = defClusters.Histo1D(
        {"hThetaAll", "#theta distribution of all tracks;#theta [rad];Counts", 100, 0, 180}, "fThetaLight");
    auto hSizeVoxelsLight = defClusters.Histo1D(
        {"hSizeVoxelsLight", "Size of light particle clusters after continuity;# Clusters;Counts", 10, 0, 10},
        "nClustersInLight");

    auto hQaveLight = defClusters.Histo1D(
        {"hQaveLight", "Average charge of light particle clusters after continuity;# Clusters;Counts", 100, 0, 10000},
        "fLight.fQave");

    // Plotthe PID of those particles to see where they lay:
    auto hPID = defClustersFiltered.Histo2D(
        {"hPID", "PID distribution of uncontinuous tracks;TL;Qtot", 200, 0, 120, 2000, 0, 3e5}, "fLight.fRawTL", "fLight.fQtotal");

    // Get some events to inspect
    std::ofstream outEvent("./Outputs/events0Cluster.dat");
    defClusters0Filtered.Foreach([&outEvent](ActRoot::MergerData& m) { m.Stream(outEvent); }, {"MergerData"});

    auto* c0 {new TCanvas {"c0", "Uncontinuous tracks"}};
    c0->Divide(2, 2);
    c0->cd(1);
    hPhi->DrawClone();
    c0->cd(2);
    hTheta->DrawClone();
    c0->cd(3);
    hPhiAll->DrawClone();
    c0->cd(4);
    hThetaAll->DrawClone();

    auto* c1 {new TCanvas {"c1", "Size of light particle clusters after continuity"}};
    hSizeVoxelsLight->DrawClone();
    auto* c2 {new TCanvas {"c2", "Average charge of light particle clusters after continuity"}};
    hQaveLight->DrawClone();
    auto* c3 {new TCanvas {"c3", "PID distribution of uncontinuous tracks"}};
    hPID->DrawClone("colz");
}