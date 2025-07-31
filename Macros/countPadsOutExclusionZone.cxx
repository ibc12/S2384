#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActSilData.h"
#include "ActTPCData.h"
#include "ActTPCParameters.h"
#include "ActCutsManager.h"

#include <ROOT/RDataFrame.hxx>

#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"

#include <fstream>
#include <iostream>

bool IsOutsideBeamZone(ROOT::Math::XYZPointF point, ActRoot::TPCParameters &tpcPars)
{
    int nPadOffset{8};
    bool isInY{point.Y() > (tpcPars.GetNPADSY() / 2 + nPadOffset) || point.Y() < (tpcPars.GetNPADSY() / 2 - nPadOffset)};
    return (isInY);
}

void countPadsOutExclusionZone()
{
    ActRoot::DataManager dataman{"../configs/data.conf", ActRoot::ModeType::EMerge};
    auto chain{dataman.GetChain()};
    // Add friends if necessary
    auto friend1{dataman.GetChain(ActRoot::ModeType::EReadTPC)};
    chain->AddFriend(friend1.get());
    auto friend2{dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(friend2.get());

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df{*chain};

    auto tpcPars{ActRoot::TPCParameters("Actar")};
    // Count the number of events
    auto dfFilterL1 = df.Filter(
        [&tpcPars](ActRoot::MergerData &d, ActRoot::ModularData &mod)
        {
            if (mod.Get("GATCONF") == 8)
            {
                return true;
            }
            else
                return false;
        },
        {"MergerData", "ModularData"});

    auto dfcounts = dfFilterL1.Define("counts",
                                      [&tpcPars](ActRoot::MergerData &d, ActRoot::ModularData &mod, ActRoot::TPCData &tpc)
                                      {
                                          int count{0};
                                          for (auto &cluster : tpc.fClusters)
                                          {
                                              auto &voxels = cluster.GetRefToVoxels();
                                              for (const auto &voxel : voxels)
                                              {
                                                  if (IsOutsideBeamZone(voxel.GetPosition(), tpcPars))
                                                      count++;
                                              }
                                          }

                                          return count;
                                      },
                                      {"MergerData", "ModularData", "TPCData"});

    auto histCounts = dfcounts.Histo1D({"counts", "Pads outside exclusion zone;Count", 300, 0, 300}, "counts");
    auto canvas{new TCanvas("canvas", "Pads Outside Exclusion Zone")};
    histCounts->DrawClone();
}