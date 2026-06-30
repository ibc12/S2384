#include "ActContinuity.h"
#include "ActMergerData.h"
#include "ActTPCData.h"
#include "ActTPCParameters.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"

void BrokenTracksZL1()
{
    // ROOT::EnableImplicitMT();

    ROOT::RDataFrame df {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_7Li_d_d_filtered.root"};
    auto def {
        df.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"})}; // only L1

    // Print the number of entries in the dataframe
    auto nEntries = static_cast<int>(def.Count().GetValue());
    std::cout << "Number of entries in the dataframe: " << nEntries << std::endl;

    // Get tpc pars
    ActRoot::TPCParameters tpc {"Actar"};

    // Build continuity object
    auto cont = std::make_shared<ActAlgorithm::Continuity>(&tpc, 20);
    // Go through the light particle and filter uncontinuous tracks
    auto defClusters = def.Define("nClustersInLight",
               [&](ActRoot::TPCData& tpc, ActRoot::MergerData& m)
               {
                   auto lightIdx = m.fLightIdx;
                   auto lightCl = tpc.fClusters[lightIdx];
                   auto voxels = lightCl.GetVoxels();
                   // Run continuity algorithm
                   auto results = cont->Run(voxels);
                   return results.first.size();
               },
               {"TPCData", "MergerData"});
    def.Foreach(
        [&](ActRoot::TPCData& tpc, ActRoot::MergerData& m)
        {
            auto lightIdx = m.fLightIdx;
            auto lightCl = tpc.fClusters[lightIdx];
            auto voxels = lightCl.GetVoxels();
            // Run continuity algorithm
            auto results = cont->Run(voxels);
            if(results.first.size() > 1)
                return true; // uncontinuous track
            else
                return false; // continuous track
        },
        {"TPCData", "MergerData"});

    // Plot the phi and theta distributions of the uncontinuous tracks
    auto hPhi = def.Histo1D({"hPhi", "#phi distribution of uncontinuous tracks;#phi [rad];Counts", 100, -180, 180},
                            "fPhiLight");
    auto hTheta = def.Histo1D({"hTheta", "#theta distribution of uncontinuous tracks;#theta [rad];Counts", 100, 0, 180},
                              "fThetaLight");
    auto hSizeVoxelsLight = defClusters.Histo1D(
        {"hSizeVoxelsLight", "Size of light particle voxels;# voxels;Counts", 10, 0, 10}, "nClustersInLight");

    auto* c0 {new TCanvas {"c0", "Uncontinuous tracks"}};
    c0->Divide(3, 1);
    c0->cd(1);
    hPhi->DrawClone();
    c0->cd(2);
    hTheta->DrawClone();
    c0->cd(3);
    hSizeVoxelsLight->DrawClone();
}