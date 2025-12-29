#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"

#include <fstream>

#include "../PrettyStyle.C"

void GateEventsFrontSiliconsVerticalShift()
{
    PrettyStyle(true);
    std::string dataconf {};
    dataconf = "./../configs/data.conf";

    // Read data
    ActRoot::DataManager dataman {dataconf, ActRoot::ModeType::EMerge};
    dataman.SetRuns(68, 122); // L0 trigger window changed for
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    auto chain3 {dataman.GetChain(ActRoot::ModeType::EFilter)};
    auto chain4 {dataman.GetChain(ActRoot::ModeType::EReadTPC)};
    chain->AddFriend(chain2.get());
    chain->AddFriend(chain3.get(), "TPCData");
    chain->AddFriend(chain4.get(), "GETTree");

    ROOT::RDataFrame df {*chain};

    // Silicons events to front
    auto dfFilterBad =
        df.Filter([](ActRoot::MergerData& m, ActRoot::ModularData& mod) { return (m.fLight.IsFilled() == true); },
                  {"MergerData", "ModularData"})
            .Filter(
                [](ActRoot::MergerData& m)
                {
                    // Gate on sil index == 2 and with z of sp > 270
                    return (m.fLight.fLayers.front() == "f0" && m.fLight.fNs.front() == 2 && m.fLight.fSP.Z() > 270.0);
                },
                {"MergerData"});

    auto dfFilterGood =
        df.Filter([](ActRoot::MergerData& m, ActRoot::ModularData& mod) { return (m.fLight.IsFilled() == true); },
                  {"MergerData", "ModularData"})
            .Filter(
                [](ActRoot::MergerData& m)
                {
                    // Gate on sil index == 2 and with z of sp > 270
                    return (m.fLight.fLayers.front() == "f0" && m.fLight.fNs.front() == 2 && m.fLight.fSP.Z() < 230.0);
                },
                {"MergerData"});

    // Forech to the events with output
    std::ofstream out("./Outputs/eventsWithShiftOnFrontLayerSil2Bad.dat");
    dfFilterBad.Foreach([&](ActRoot::MergerData& m) { m.Stream(out); }, {"MergerData"});
    out.close();

    std::ofstream out1("./Outputs/eventsWithShiftOnFrontLayerSil2Good.dat");
    dfFilterGood.Foreach([&](ActRoot::MergerData& m) { m.Stream(out1); }, {"MergerData"});
    out1.close();
}