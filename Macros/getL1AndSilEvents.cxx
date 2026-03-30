#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"

#include "ROOT/RDataFrame.hxx"

#include <fstream>
#include <string>


void getL1AndSilEvents()
{
    // Get data 
    ActRoot::DataManager dataman {"../configs/data_7Li_auxiliar.conf", ActRoot::ModeType::EMerge};
    // Select just 4 runs, 2 of 11Li, and 2 of 7Li
    auto chain {dataman.GetChain()};
    auto chain1 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(chain1.get());
    ROOT::RDataFrame df {*chain};

    auto dfL1 = df.Filter(
        [](ActRoot::ModularData& m)
        {
            if(m.Get("GATCONF") == 8)
                return true;
            return false;
        },
        {"ModularData"});
    auto dfL1WithLight = df.Filter(
        [](ActRoot::ModularData& m, ActRoot::MergerData& mer)
        {
            if(m.Get("GATCONF") == 8 && mer.fLightIdx != -1)
                return true;
            return false;
        },
        {"ModularData", "MergerData"});
    auto dfSil = df.Filter(
        [](ActRoot::ModularData& m)
        {
            if(m.Get("GATCONF") == 1 || m.Get("GATCONF") == 2)
                return true;
            return false;
        },
        {"ModularData"});

    auto dfSilFront = df.Filter(
        [](ActRoot::ModularData& m)
        {
            if(m.Get("GATCONF") == 4)
                return true;
            return false;
        },
        {"ModularData"});

    // Get events into file .dat
    // std::ofstream streamer {"./Outputs/L1_run66_to_82.dat"};
    // dfL1.Foreach(
    //     [&](ActRoot::MergerData& m)
    //     {
    //         m.Stream(streamer);
    //     },
    //     {"MergerData"});
    // streamer.close();
    std::ofstream streamer1 {"./Outputs/L1WithLight_run66_to_82.dat"};
    dfL1WithLight.Foreach(
        [&](ActRoot::MergerData& m)
        {
            m.Stream(streamer1);
        },
        {"MergerData"});
    streamer1.close();
    // std::ofstream streamer2 {"./Outputs/Sil_run66_to_82.dat"};
    // dfSil.Foreach(
    //     [&](ActRoot::MergerData& m)
    //     {
    //         m.Stream(streamer2);
    //     },
    //     {"MergerData"});
    // streamer2.close();
    // std::ofstream streamer3 {"./Outputs/SilFront_run66_to_82.dat"};
    // dfSilFront.Foreach(
    //     [&](ActRoot::MergerData& m)
    //     {
    //         m.Stream(streamer3);
    //     },
    //     {"MergerData"});
    // streamer3.close();

}