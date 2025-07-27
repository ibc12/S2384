#include "ActDataManager.h"
#include "ActSilData.h"
#include "ActTPCData.h"
#include "ActModularData.h"
#include "ActModularData.h"
#include "ActMergerData.h"
#include "ActTPCData.h"

#include <ROOT/RDataFrame.hxx>

#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"

void BeamCountsOverTime()
{
    // Get Sil data
    ActRoot::DataManager dataman{"../configs/data.conf", ActRoot::ModeType::EReadSilMod};
    auto chain{dataman.GetChain()};

    auto friend1{dataman.GetChain(ActRoot::ModeType::EFilter)};
    chain->AddFriend(friend1.get());

    // Build the RDataFrame
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df{*chain};

    // Gate on beam events, in principle CFA trigger
    auto gated{df.Filter([](ActRoot::ModularData &d)
                         { return d.Get("GATCONF") == 64; }, {"ModularData"})
                   .Filter(
                       [](ActRoot::TPCData &d)
                       {
                           return d.fClusters.size() >= 1;
                       },
                       {"TPCData"})
                   .Filter([](ActRoot::TPCData &d)
                           {
            auto [xmin, xmax] = d.fClusters.front().GetXRange();
            return (xmin < 5 && xmax > 120); }, {"TPCData"})};

    // Define the size variable
    gated = gated.Define("size", [](const ActRoot::TPCData &d)
                         { return std::count_if(d.fClusters.begin(), d.fClusters.end(), [](const auto &clu)
                                                {
            auto [xmin, xmax] = clu.GetXRange();
            return (xmin < 5 && xmax > 120); }); }, {"TPCData"});

    auto histo = gated.Histo1D({"hBeamCounts", "Beam Counts DIstribution;Number of beam particles in an event;Counts", 7, 0, 7}, "size");
    // Create canvas and draw histogram
    TCanvas *c1 = new TCanvas("c1", "Beam Counts Over Time", 800, 600);
    histo->DrawClone();
}