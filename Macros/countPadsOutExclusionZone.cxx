#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActSilData.h"
#include "ActTPCData.h"
#include "ActTPCParameters.h"

#include <ROOT/RDataFrame.hxx>

#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"

#include <fstream>
#include <iostream>

bool IsOutsideBeamZone(ROOT::Math::XYZPointF point, ActRoot::TPCParameters& tpcPars)
{
    int nPadMin {56};
    int nPadMax {71};
    bool isInY {point.Y() > (nPadMax) || point.Y() < (nPadMin)};
    return (isInY);
}

void countPadsOutExclusionZone()
{
    ActRoot::DataManager dataman {"../configs/data.conf", ActRoot::ModeType::EMerge};
    auto chain {dataman.GetChain()};
    // Add friends if necessary
    auto friend1 {dataman.GetChain(ActRoot::ModeType::EReadTPC)};
    chain->AddFriend(friend1.get());
    auto friend2 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(friend2.get());

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {*chain};

    auto tpcPars {ActRoot::TPCParameters("Actar")};
    // Count the number of events
    auto dfFilterL1 = df.Filter(
        [&tpcPars](ActRoot::MergerData& mer, ActRoot::ModularData& mod)
        {
            if(mod.Get("GATCONF") == 8 && (mer.fLightIdx != -1))
            {
                return true;
            }
            else
                return false;
        },
        {"MergerData", "ModularData"});

    auto dfcounts = dfFilterL1.Define("nPadsOutside",
                                      [&tpcPars](ActRoot::TPCData& tpc)
                                      {
                                          std::set<std::pair<int, int>> pads;

                                          for(auto& cluster : tpc.fClusters)
                                          {
                                              auto& voxels = cluster.GetRefToVoxels();
                                              for(const auto& voxel : voxels)
                                              {
                                                  const auto& pos = voxel.GetPosition();

                                                  if(IsOutsideBeamZone(pos, tpcPars))
                                                  {
                                                      int padx = static_cast<int>(std::floor(pos.X()));
                                                      int pady = static_cast<int>(std::floor(pos.Y()));
                                                      pads.insert({padx, pady});
                                                  }
                                              }
                                          }
                                          for(auto& voxels : tpc.fRaw)
                                          {
                                              const auto& pos = voxels.GetPosition();

                                              if(IsOutsideBeamZone(pos, tpcPars))
                                              {
                                                  int padx = static_cast<int>(std::floor(pos.X()));
                                                  int pady = static_cast<int>(std::floor(pos.Y()));
                                                  pads.insert({padx, pady});
                                              }
                                          }

                                          return static_cast<int>(pads.size());
                                      },
                                      {"TPCData"});
    auto dfBeforeL1Thresholdchange =
        dfcounts.Filter([](ActRoot::MergerData& m) { return m.fRun < 58; }, {"MergerData"});
    auto dfAfterL1Thresholdchange = dfcounts.Filter([](ActRoot::MergerData& m) { return m.fRun > 58; }, {"MergerData"});
    auto histCounts =
        dfcounts.Histo1D({"nPadsOutside", "Pads outside exclusion zone;Count", 300, 0, 300}, "nPadsOutside");
    auto histCountsBeforeL1Thresholdchange =
        dfBeforeL1Thresholdchange.Histo1D({"nPadsOutsideBeforeL1Thresholdchange",
                                           "Pads outside exclusion zone before L1 threshold change;Count", 300, 0, 300},
                                          "nPadsOutside");
    auto histCountsAfterL1Thresholdchange =
        dfAfterL1Thresholdchange.Histo1D({"nPadsOutsideAfterL1Thresholdchange",
                                          "Pads outside exclusion zone after L1 threshold change;Count", 300, 0, 300},
                                         "nPadsOutside");
    auto canvas {new TCanvas("canvas", "Pads Outside Exclusion Zone")};
    canvas->Divide(3, 1);
    canvas->cd(1);
    histCounts->DrawClone();
    canvas->cd(2);
    histCountsBeforeL1Thresholdchange->DrawClone();
    canvas->cd(3);
    histCountsAfterL1Thresholdchange->DrawClone();

    // Save events if needed
    // std::ofstream out(TString::Format("./Outputs/eventsLessThan4PadsOutExclusion.dat").Data());
    // dfcounts.Foreach(
    //     [&](ActRoot::MergerData& m, int nPadsOutside)
    //     {
    //         if(nPadsOutside < 4)
    //             m.Stream(out);
    //     },
    //     {"MergerData", "nPadsOutside"});
    // out.close();
}