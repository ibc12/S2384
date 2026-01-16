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
    int nPadOffset {8};
    bool isInY {point.Y() > (tpcPars.GetNPADSY() / 2 + nPadOffset) ||
                point.Y() < (tpcPars.GetNPADSY() / 2 - nPadOffset)};
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

    //ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {*chain};

    auto tpcPars {ActRoot::TPCParameters("Actar")};
    // Count the number of events
    auto dfFilterL1 = df.Filter(
        [&tpcPars](ActRoot::MergerData& d, ActRoot::ModularData& mod)
        {
            if(mod.Get("GATCONF") == 8)
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
                                                      int padx = static_cast<int>(std::floor(pos.X() / 2.0));
                                                      int pady = static_cast<int>(std::floor(pos.Y() / 2.0));
                                                      pads.insert({padx, pady});
                                                  }
                                              }
                                          }

                                          return static_cast<int>(pads.size());
                                      },
                                      {"TPCData"});
    auto histCounts =
        dfcounts.Histo1D({"nPadsOutside", "Pads outside exclusion zone;Count", 300, 0, 300}, "nPadsOutside");
    auto canvas {new TCanvas("canvas", "Pads Outside Exclusion Zone")};
    histCounts->DrawClone();

    // Save events if needed
    std::ofstream out(TString::Format("./Outputs/eventsOnePadOutExclusion.dat").Data());
    dfcounts.Foreach(
        [&](ActRoot::MergerData& m, int nPadsOutside)
        {
            if(nPadsOutside > 50)
                m.Stream(out);
        },
        {"MergerData", "nPadsOutside"});
    out.close();
}