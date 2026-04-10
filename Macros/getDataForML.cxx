#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActParticle.h"
#include "ActSRIM.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TH1.h"
#include "TMath.h"
#include "TString.h"

#include <fstream>
#include <iostream>
#include <string>


void getDataForML()
{
    std::string infile {"../PostAnalysis/Outputs/tree_preprocess_F_7Li.root"};
    std::string infile1 {"../PostAnalysis/Outputs/tree_preprocess_7Li.root"};

    // RDataFrame
    // ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"PreProcessed_Tree", infile};
    ROOT::RDataFrame df1 {"PreProcessed_Tree", infile1};

    // Get drift parameter from config file
    ActRoot::InputParser parser {};
    parser.ReadFile("../configs/detector.conf");
    auto driftBlock = parser.GetBlock("Merger");
    auto driftFactor = driftBlock->GetDouble("DriftFactor"); // in mm^2/us

    auto dfFiltered = df.Filter([](ActRoot::MergerData& mer, ActRoot::ModularData& mod)
              { return (mod.Get("GATCONF") == 8) && (mer.fLightIdx != -1); }, {"MergerData", "ModularData"});
    auto df1Filtered = df1.Filter([](ActRoot::MergerData& mer, ActRoot::ModularData& mod)
               { return (mod.Get("GATCONF") == 8) && (mer.fLightIdx != -1); }, {"MergerData", "ModularData"});

    // Go through all events, and save the information for the voxels in a root file
    auto dfDefines =
        dfFiltered.Define("X",
                  [](ActRoot::MergerData& mer, ActRoot::TPCData& tpc)
                  {
                      ROOT::RVec<float> x;
                      if(mer.fLightIdx < 0 || static_cast<std::size_t>(mer.fLightIdx) >= tpc.fClusters.size())
                          return x;
                      auto& cluster = tpc.fClusters.at(mer.fLightIdx);
                      auto& voxels = cluster.GetRefToVoxels();

                      for(const auto& v : voxels)
                      {
                          auto pos = v.GetPosition();
                          x.push_back(pos.X() * 2);
                      }
                      return x;
                  },
                  {"MergerData", "TPCData"})
            .Define("Y",
                    [](ActRoot::MergerData& mer, ActRoot::TPCData& tpc)
                    {
                        ROOT::RVec<float> y;
                        if(mer.fLightIdx < 0 || static_cast<std::size_t>(mer.fLightIdx) >= tpc.fClusters.size())
                            return y;
                        auto& cluster = tpc.fClusters.at(mer.fLightIdx);
                        auto& voxels = cluster.GetRefToVoxels();

                        for(const auto& v : voxels)
                        {
                            auto pos = v.GetPosition();
                            y.push_back(pos.Y() * 2);
                        }
                        return y;
                    },
                    {"MergerData", "TPCData"})
            .Define("Z",
                    [&driftFactor](ActRoot::MergerData& mer, ActRoot::TPCData& tpc)
                    {
                        ROOT::RVec<float> z;
                        if(mer.fLightIdx < 0 || static_cast<std::size_t>(mer.fLightIdx) >= tpc.fClusters.size())
                            return z;
                        auto& cluster = tpc.fClusters.at(mer.fLightIdx);
                        auto& voxels = cluster.GetRefToVoxels();

                        for(const auto& v : voxels)
                        {
                            auto pos = v.GetPosition();
                            z.push_back(pos.Z() * driftFactor);
                        }
                        return z;
                    },
                    {"MergerData", "TPCData"})
            .Define("Charge",
                    [](ActRoot::MergerData& mer, ActRoot::TPCData& tpc)
                    {
                        ROOT::RVec<float> charge;
                        if(mer.fLightIdx < 0 || static_cast<std::size_t>(mer.fLightIdx) >= tpc.fClusters.size())
                            return charge;
                        auto& cluster = tpc.fClusters.at(mer.fLightIdx);
                        auto& voxels = cluster.GetRefToVoxels();

                        for(const auto& v : voxels)
                            charge.push_back(v.GetCharge());
                        return charge;
                    },
                    {"MergerData", "TPCData"});

    dfDefines.Snapshot("ML_Tree", "./Outputs/dataForML.root", {"X", "Y", "Z", "Charge"});

    int counter = 0;
    // Print the voxel positions and charge for all the events
    dfDefines.Foreach(
        [&counter](ROOT::RVec<float> x, ROOT::RVec<float> y, ROOT::RVec<float> z, ROOT::RVec<float> charge)
        {
            if(counter > 10)
            {
            }
            else
            {
                std::cout << "========= EVENT =========\n";
                for(size_t i = 0; i < x.size(); ++i)
                {
                    std::cout << "Voxel " << i << ": X = " << x[i] << " mm, Y = " << y[i] << " mm, Z = " << z[i]
                              << " mm, Charge = " << charge[i] << "\n";
                }
                counter++;
            }
        },
        {"X", "Y", "Z", "Charge"});

    auto hZdrift = dfFiltered.Histo1D({"hZdrift", "Z drift distribution;Z drift [mm];Entries", 350, -50, 300}, "zDrift");
    auto hZ = dfDefines.Histo1D({"hZ", "Z distribution;Z [mm];Entries", 500, -50, 450}, "Z");

    auto hZdrift1 = df1Filtered.Histo1D({"hZdrift1", "Z drift distribution;Z drift [mm];Entries", 350, -50, 300}, "zDrift");

    auto* c1 = new TCanvas("c1", "c1", 1200, 600);
    c1->Divide(3, 1);
    c1->cd(1);
    hZ->DrawClone();
    c1->cd(2);
    hZdrift->DrawClone();
    c1->cd(3);
    hZdrift1->DrawClone();
}