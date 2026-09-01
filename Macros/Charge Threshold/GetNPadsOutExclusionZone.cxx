#include "ActMergerData.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH2D.h>
#include <TPaveText.h>
#include <TProfile.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

constexpr int yMinExclusionZone = 56;
constexpr int yMaxExclusionZone = 71;

void GetNPadsOutExclusionZone()
{
    ROOT::EnableImplicitMT();
    std::string file {"../../PostAnalysis/Outputs/tree_preprocess_F_7Li.root"};

    ROOT::RDataFrame df_exp {"PreProcessed_Tree", file.c_str()};

    // Filter to just end with r66e292 r66e973 r66e1554 r66e6583 r66e8012 r67e23132
    auto df_exp_filtered = df_exp.Filter(
        [](ActRoot::MergerData& m)
        {
            return (m.fRun == 66 && m.fEntry == 292) || (m.fRun == 66 && m.fEntry == 973) ||
                   (m.fRun == 66 && m.fEntry == 1554) || (m.fRun == 66 && m.fEntry == 6583) ||
                   (m.fRun == 66 && m.fEntry == 8012) || (m.fRun == 67 && m.fEntry == 23132);
        },
        {"MergerData"});

    // Define nPads outside exclusion zone
    auto df_exp_nPads = df_exp_filtered.Define("nPads",
                                               [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                                               {
                                                   int nPads {};
                                                   for(auto cluster : tpc.fClusters)
                                                   {
                                                       for(auto voxel : cluster.GetRefToVoxels())
                                                       {
                                                           int iy = voxel.GetPosition().Y();
                                                           if(iy < yMinExclusionZone || iy > yMaxExclusionZone)
                                                               nPads++;
                                                       }
                                                   }
                                                   return nPads;
                                               },
                                               {"MergerData", "TPCData"});

    // Print the the run and entry number together with the nPads value for each event
    df_exp_nPads.Foreach(
        [](ActRoot::MergerData& m, int nPads)
        { std::cout << "Run: " << m.fRun << ", Entry: " << m.fEntry << ", nPads: " << nPads << std::endl; },
        {"MergerData", "nPads"});
}