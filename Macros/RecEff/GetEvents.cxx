#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"

#include "ROOT/RDataFrame.hxx"

#include <fstream>

void GetEvents()
{
    // Get the runs from the data manager reading the config file
    ActRoot::DataManager dataManager {};
    dataManager.ReadDataFile("../../configs/data.conf");

    // For L1 enough with 4 runs (64, 67). For lat sils put at least 10
    dataManager.SetRuns(50, 51);

    // Get df for the runs
    auto chain {dataManager.GetChain(ActRoot::ModeType::EReadSilMod)};
    auto chain1 {dataManager.GetChain(ActRoot::ModeType::EMerge)};
    chain->AddFriend(chain1.get());
    auto df {ROOT::RDataFrame(*chain)};

    // GATCONFs: 1 -> l0, 2 -> r0, 4 -> f0, 8 -> L1,

    auto dfFilterGATCONF = df.Filter(
        [](ActRoot::ModularData& mod)
        {
            if(mod.Get("GATCONF") == 8)
                return true;
            return false;
        },
        {"ModularData"});

    std::ofstream outFile("./Inputs/events_L1_11Li_preChange.dat");
    // Save events
    dfFilterGATCONF.Foreach([&outFile](ActRoot::MergerData& m) { m.Stream(outFile); }, {"MergerData"});
}