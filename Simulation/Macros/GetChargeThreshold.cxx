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
#include "TString.h"

#include <fstream>
#include <iostream>
#include <map>
#include <string>

#include "./../../PrettyStyle.C"

void GetChargeThreshold()
{
    // Get some runs of data and analyze the charge threshold
    // It has to be done without pad align
    // PrettyStyle(false);
    std::string dataconf {"../../configs/data.conf"};

    // Read data
    ActRoot::DataManager dataman {dataconf, ActRoot::ModeType::EReadTPC};
    dataman.SetRuns(64, 67);
    auto chain {dataman.GetChain()};
    // auto chain2 {dataman.GetChain(ActRoot::ModeType::EMerge)};
    // chain->AddFriend(chain2.get());

    // RDataFrame
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {*chain};

    auto dfFilter = df.Filter(
        [](ActRoot::TPCData& tpc)
        {
            if(tpc.fClusters.size() == 0 && tpc.fRaw.size() == 0)
                return false;
            else
                return true;
        },
        {"TPCData"});

    auto dfThresh = dfFilter.Define("chargeThreshold",
                                    [](ActRoot::TPCData& tpc)
                                    {
                                        float minCharge = 1e6;
                                        for(auto cluster : tpc.fClusters)
                                        {
                                            for(auto voxel : cluster.GetVoxels())
                                            {
                                                float charge = voxel.GetCharge();
                                                if(charge < minCharge)
                                                    minCharge = charge;
                                            }
                                        }

                                        for(auto voxel : tpc.fRaw)
                                        {
                                            float charge = voxel.GetCharge();
                                            if(charge < minCharge)
                                                minCharge = charge;
                                        }

                                        return minCharge;
                                    },
                                    {"TPCData"});

    auto dfMax = dfFilter.Define("maxCharge",
                                 [](ActRoot::TPCData& tpc)
                                 {
                                     float maxCharge = -10;
                                     for(auto cluster : tpc.fClusters)
                                     {
                                         for(auto voxel : cluster.GetVoxels())
                                         {
                                             float charge = voxel.GetCharge();
                                             if(charge > maxCharge && voxel.GetIsSaturated() == true)
                                                 maxCharge = charge;
                                         }
                                     }

                                     for(auto voxel : tpc.fRaw)
                                     {
                                         float charge = voxel.GetCharge();
                                         if(charge > maxCharge && voxel.GetIsSaturated() == true)
                                             maxCharge = charge;
                                     }

                                     return maxCharge;
                                 },
                                 {"TPCData"});

    auto hThresh = dfThresh.Histo1D({"chargeThreshold", "Charge Threshold", 1000, 0, 1000}, "chargeThreshold");
    hThresh->GetXaxis()->SetTitle("Charge Threshold [channels]");
    hThresh->GetYaxis()->SetTitle("Entries");
    auto hMax = dfMax.Histo1D({"maxCharge", "Max Charge", 10000, 0, 10000}, "maxCharge");
    hMax->GetXaxis()->SetTitle("Max Charge [channels]");
    hMax->GetYaxis()->SetTitle("Entries");

    // Save events if needed
    // std::ofstream out(TString::Format("./Outputs/eventsNegativeMinCharge.dat").Data());
    // dfThresh.Foreach(
    //     [&](ActRoot::MergerData& m, float charge)
    //     {
    //         if(charge < 0)
    //             m.Stream(out);
    //     },
    //     {"MergerData", "chargeThreshold"});
    // out.close();
    // std::ofstream out1(TString::Format("./Outputs/eventsMaxCharge.dat").Data());
    // dfMax.Foreach(
    //     [&](ActRoot::MergerData& m, float charge)
    //     {
    //         if(charge > 7000)
    //             m.Stream(out1);
    //     },
    //     {"MergerData", "maxCharge"});
    // out1.close();

    auto* c1 = new TCanvas();
    c1->Divide(2, 1);
    c1->cd(1);
    hThresh->DrawClone();
    c1->cd(2);
    hMax->DrawClone();
}