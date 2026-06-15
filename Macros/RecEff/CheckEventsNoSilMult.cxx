#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActSilData.h"

#include "ROOT/RDataFrame.hxx"

#include <fstream>

void CheckEventsNoSilMult()
{
    // Get the runs from the data manager reading the config file
    ActRoot::DataManager dataManager {};
    dataManager.ReadDataFile("../../configs/data.conf");

    // For L1 enough with 4 runs (64, 67). For lat sils put at least 10
    dataManager.SetRuns(60, 70);

    // Get df for the runs
    auto chain {dataManager.GetChain(ActRoot::ModeType::EReadSilMod)};
    auto chain1 {dataManager.GetChain(ActRoot::ModeType::EMerge)};
    chain->AddFriend(chain1.get());
    auto df {ROOT::RDataFrame(*chain)};

    auto dfEvent = df.Filter([](ActRoot::MergerData& m) { return m.fRun == 60 && m.fEntry == 21626; }, {"MergerData"});

    dfEvent.Foreach([](ActRoot::MergerData& m) { m.Print(); }, {"MergerData"});

    dfEvent.Foreach(
        [](ActRoot::SilData& sil)
        {
            auto layers = sil.GetLayers();
            for(auto& layer : layers)
            {
                std::cout << "Layer: " << layer << " with mult: " << sil.GetMult(layer) << std::endl;
                for(int i = 0; i < sil.GetMult(layer); i++)
                {
                    auto detID = sil.fSiN[layer][i];
                    auto energy = sil.fSiE[layer][i];
                    std::cout << "  DetID: " << detID << ", Energy: " << energy << " keV" << std::endl;
                }
            }
        },
        {"SilData"});
}