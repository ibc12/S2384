#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActSilData.h"

#include "ROOT/RDataFrame.hxx"

#include <fstream>

void DebugF2Multiplicity()
{
    // This Macro do no work in this directory, but a directory up it does ?¿?¿?¿



    // Get the runs from the data manager reading the config file
    ActRoot::DataManager dataManager {};
    dataManager.ReadDataFile("../../../configs/data.conf");

    // For L1 enough with 4 runs (64, 67). For lat sils put at least 10
    dataManager.SetRuns(60, 70);

    // Get df for the runs
    auto chain {dataManager.GetChain(ActRoot::ModeType::EReadSilMod)};
    auto chain1 {dataManager.GetChain(ActRoot::ModeType::EMerge)};
    chain->AddFriend(chain1.get(), ("MergerData"));
    auto df {ROOT::RDataFrame(*chain)};

    // GATCONFs: 1 -> l0, 2 -> r0, 4 -> f0, 8 -> L1,

    df.Describe().Print();

    // auto dfFilterGATCONF = df.Filter(
    //     [](ActRoot::ModularData& mod)
    //     {
    //         if(mod.Get("GATCONF") == 1 || mod.Get("GATCONF") == 2)
    //             return true;
    //         return false;
    //     },
    //     {"ModularData"});
    //
    // // Definethe multiplicity of F2
    // auto dfMultiplicityF2 = dfFilterGATCONF.Define(
    //     "F2Multiplicity",
    //     [](ActRoot::SilData& s)
    //     {
    //         for(auto layer : s.GetLayers())
    //         {
    //             if(layer == "f2")
    //             {
    //                 int size = s.fSiN[layer].size();
    //                 return size;
    //             }
    //         }
    //         return -1111;
    //     },
    //     {"SilData"});
    //
    // // Get events with F2 multiplicity = 1
    // auto dfFilterF2Multiplicity1 = dfMultiplicityF2.Filter(
    //     [](ActRoot::MergerData& m, ActRoot::SilData& s)
    //     {
    //         for(auto layer : s.GetLayers())
    //         {
    //             if(layer == "f2")
    //             {
    //                 if(s.fSiN[layer].size() == 1)
    //                     return true;
    //                 else
    //                     return false;
    //             }
    //         }
    //         return false;
    //     },
    //     {"MergerData", "SilData"});
    // auto dfFitlerF2MultiplicityGreater1 = dfMultiplicityF2.Filter(
    //     [](ActRoot::MergerData& m, ActRoot::SilData& s)
    //     {
    //         for(auto layer : s.GetLayers())
    //         {
    //             if(layer == "f2")
    //             {
    //                 if(s.fSiN[layer].size() > 1)
    //                     return true;
    //                 else
    //                     return false;
    //             }
    //         }
    //         return false;
    //     },
    //     {"MergerData", "SilData"});
    //
    // // Plot the multiplicity of F2
    // auto hMultiplicityF2 = dfMultiplicityF2.Histo1D(
    //     {"hMultiplicityF2", "Multiplicity of F2;Multiplicity;Counts", 5, -0.5, 4.5},
    //     "F2Multiplicity");
    //
    // // Plot the multiplicity of F2 for events with F2 multiplicity > 1
    // auto hMultiplicityF2Greater1 = dfFitlerF2MultiplicityGreater1.Histo1D(
    //     {"hMultiplicityF2Greater1", "Multiplicity of F2 for events with F2 multiplicity > 1;Multiplicity;Counts", 5,
    //     -0.5, 4.5}, "F2Multiplicity");
}