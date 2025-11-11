#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"

#include "ROOT/RDataFrame.hxx"

#include <fstream>
#include <string>


void getL1AndSilEvents()
{
    // Get data 
    ActRoot::DataManager dataman {"../configs/data.conf", ActRoot::ModeType::EMerge};
    // Select just 4 runs, 2 of 11Li, and 2 of 7Li
    dataman.SetRuns(62, 69);
    auto chain {dataman.GetChain()};
    auto chain1 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(chain1.get());
    ROOT::RDataFrame df {*chain};

    auto dfL1AndSil = df.Filter(
        [](ActRoot::ModularData& m)
        {
            if(m.Get("GATCONF") == 8)
                return true;
            return false;
        },
        {"ModularData"});

    // Get events into file .dat
    std::ofstream streamer {"./Outputs/L1_run62_to_69.dat"};
    dfL1AndSil.Foreach(
        [&](ActRoot::MergerData& m)
        {
            m.Stream(streamer);
        },
        {"MergerData"});
    streamer.close();

}