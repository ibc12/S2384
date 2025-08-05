#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActTPCData.h"
#include "ActTypes.h"

#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TROOT.h"

#include <fstream>
#include <iostream>
void l1_with_range()
{
    ActRoot::DataManager dataman {"../configs/data.conf", ActRoot::ModeType::EMerge};
    // dataman.SetRuns(55, 64);
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    auto chain3 {dataman.GetChain(ActRoot::ModeType::EFilter)};
    chain->AddFriend(chain2.get());
    chain->AddFriend(chain3.get());

    // Dataframe
    // ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {*chain};

    auto gated {
        df.Filter(
              [](ActRoot::MergerData& mer, ActRoot::ModularData& mod)
              {
                  if(mod.Get("GATCONF") != 8)
                      return false;
                  if(mer.fLightIdx == -1)
                      return false;
                  if(mer.fLight.HasSP())
                      return false;
                  return true;
              },
              {"MergerData", "ModularData"})
            .Define("range", [](ActRoot::MergerData& mer) { return (mer.fBraggP - mer.fRP).R(); }, {"MergerData"})
            .Define("qtotal", [](ActRoot::MergerData& mer) { return mer.fLight.fQtotal; }, {"MergerData"})
            .Define("tl",
                    [](ActRoot::MergerData& mer, ActRoot::TPCData& tpc)
                    {
                        auto& cluster {tpc.fClusters[mer.fLightIdx]};
                        cluster.SortAlongDir();
                        auto& voxels {cluster.GetVoxels()};
                        auto begin {voxels.front().GetPosition()};
                        auto end {voxels.back().GetPosition()};
                        return (begin - end).R();
                    },
                    {"MergerData", "TPCData"})};

    // Book histograms
    auto hL1PID {
        gated.Histo2D({"hL1PID", "L1 PID;Range [mm];Qtotal [au]", 240, 0, 200, 1000, 0, 300000}, "range", "qtotal")};
    auto hL1TL {gated.Histo2D({"hL1TL", "L1 PID;TL [mm];Qtotal [au]", 240, 0, 200, 1000, 0, 300000}, "tl", "qtotal")};
    auto hL1Theta {gated.Histo2D({"hL1Theta", "L1 #theta;TL [mm];#theta_{Lab} [#deg]", 240, 0, 200, 300, 0, 100}, "tl",
                                 "fThetaLight")};
    // Count
    auto count {gated.Count()};
    std::ofstream streamer {"./Outputs/L1_events.dat"};
    gated.Foreach([&](ActRoot::MergerData& mer) { mer.Stream(streamer); }, {"MergerData"});
    streamer.close();
    std::cout << "Gated : " << *count << '\n';

    // Plot
    auto* c0 {new TCanvas {"c0", "l1 plot"}};
    c0->DivideSquare(4);
    c0->cd(1);
    hL1PID->DrawClone("colz");
    c0->cd(2);
    hL1TL->DrawClone("colz");
    c0->cd(3);
    hL1Theta->DrawClone("colz");
}
