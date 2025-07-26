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

    // Stream entry number
    std::ofstream streamer {"./Outputs/gatconf_f0.dat"};
    df.Foreach(
        [&](ActRoot::ModularData& mod, ActRoot::MergerData& mer)
        {
            if(mod.Get("GATCONF") == 4) // f0
                mer.Stream(streamer);
        },
        {"ModularData", "MergerData"});
    streamer.close();
}
