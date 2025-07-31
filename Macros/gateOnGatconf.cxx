#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"

#include <fstream>

void gateOnGatconf()
{
    ActRoot::DataManager dataman {"../configs/data.conf", ActRoot::ModeType::EReadSilMod};
    auto chain {dataman.GetChain()};
    auto chainMerger {dataman.GetChain(ActRoot::ModeType::EMerge)};
    chain->AddFriend(chainMerger.get());

    ROOT::RDataFrame df {*chain};

    auto dfFilter {df.Filter(
        [](ActRoot::MergerData& m)
        {
            if(!m.fLight.IsFilled() || m.fLight.fLayers.front() != "f0")
                return false;
            if(m.fLight.fNs.front() == 2)
                return true;
            return false;
        },
        {"MergerData"})};

    // Stream entry number
    std::ofstream streamer {"./Outputs/gatconf_f0_sil2.dat"};
    dfFilter.Foreach(
        [&](ActRoot::ModularData& mod, ActRoot::MergerData& mer)
        {
            if(mod.Get("GATCONF") == 4)
                mer.Stream(streamer);
        },
        {"ModularData", "MergerData"});
    streamer.close();
}
