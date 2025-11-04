#include "ActMergerData.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TH1D.h"
#include "TString.h"

#include <fstream>
#include <iostream>
#include <string>

void DebugHighChargeAfterRP()
{
    std::string beam {"11Li"};
    std::string target {"d"};
    std::string light {"p"};

    // Get file output from pipe2
    TString filename {TString::Format("../PostAnalysis/Outputs/tree_ex_%s_%s_%s_filtered.root", beam.c_str(),
                                      target.c_str(), light.c_str())};
    ROOT::RDataFrame df {"Final_Tree", filename};

    // Get just 1 event and see the charge, for run == 20 entry 41251
    auto dfFiltered {df.Filter(
                           [](ActRoot::MergerData& m)
                           {
                               if(m.fRun == 67 && m.fEntry == 9239)
                                   return true;
                               else
                                   return false;
                           },
                           {"MergerData"})};
    
    dfFiltered.Foreach(
        [&](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
        {
            std::cout << "Run: " << m.fRun << " Entry: " << m.fEntry << "\n";
            std::cout << "RP: (" << m.fRP.X() << ", " << m.fRP.Y() << ", " << m.fRP.Z() << ")\n";
            std::cout << "Clusters:\n";
            double qMax = -1;
            for(auto& cluster : tpc.fClusters)
            {
                std::cout << " Cluster with " << cluster.GetRefToVoxels().size() << " voxels:\n";
                for(auto& voxel : cluster.GetRefToVoxels())
                {
                    auto pos {voxel.GetPosition()};
                    std::cout << "  Voxel at (" << pos.X() << ", " << pos.Y() << ", " << pos.Z() << ")"
                              << " Charge: " << voxel.GetCharge() << "\n";
                    if(voxel.GetCharge() > qMax)
                        qMax = voxel.GetCharge();
                }
            }
            std::cout << "Max charge in event: " << qMax << "\n";
        },
        {"MergerData", "GETTree_TPCData"});
}