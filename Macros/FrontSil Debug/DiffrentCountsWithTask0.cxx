#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TH2.h"
#include "TH2D.h"
#include "TString.h"

#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <mutex>
#include <unordered_set>
#include <cstdint>

void DiffrentCountsWithTask0()
{
    // Coger el run 31 de la carpeta Rootfiles/Merger y Merger_all_preTask0
    // Convertirlos a df
    // Guardar los eventos diferentes ok

    std::string inputFile = "../../RootFiles/Merger/Merged_Run_0031.root";
    std::string inputFile2 = "../../RootFiles/Merger_all_preTask0/Merged_Run_0031.root";
    std::string treeName = "ACTAR_Merged";

    std::cout << "Opening: " << inputFile << " tree: " << treeName << std::endl;
    ROOT::RDataFrame df(treeName, inputFile);

    std::cout << "Opening: " << inputFile2 << " tree: " << treeName << std::endl;
    ROOT::RDataFrame df2(treeName, inputFile2);

    auto total = *df.Count();
    std::cout << "Total entries: " << total << std::endl;

    auto total2 = *df2.Count();
    std::cout << "Total entries: " << total2 << std::endl;

    // Filter: mergerdata.fFlag == "ok"
    // Use a lambda with auto parameter so it accepts either const char* or std::string branches.
    auto df_ok = df.Filter(
        [](ActRoot::MergerData& m)
        {
            return m.fFlag.c_str() == std::string("ok");
        },
        {"MergerData"});

    auto df2_ok = df2.Filter(
        [](ActRoot::MergerData& m)
        {
            return m.fFlag.c_str() == std::string("ok");
        },
        {"MergerData"});

    
    // Get the counts after filtering
    auto count_ok = *df_ok.Count();
    std::cout << "Entries with fFlag == 'ok': " << count_ok << std::endl;  

    auto count2_ok = *df2_ok.Count();
    std::cout << "Entries with fFlag == 'ok' in second file: " << count2_ok << std::endl;

    // Build a set of (run,entry) keys from df2_ok to allow fast membership tests.
    // We'll collect keys with a thread-safe push into a vector, then build an unordered_set.
    std::vector<uint64_t> keys2;
    keys2.reserve((size_t)count2_ok);
    std::mutex keys2_mtx;

    // Foreach will materialize the values into keys2
    df2_ok.Foreach(
        [&keys2, &keys2_mtx](const ActRoot::MergerData& m)
        {
            uint64_t key = (uint64_t(m.fRun) << 32) | uint64_t(uint32_t(m.fEntry));
            std::lock_guard<std::mutex> lk(keys2_mtx);
            keys2.push_back(key);
        },
        {"MergerData"});

    // Build unordered_set for O(1) lookup
    std::unordered_set<uint64_t> keys2_set;
    keys2_set.reserve(keys2.size() * 2 + 10);
    for(auto k : keys2) keys2_set.insert(k);

    // Now filter df_ok keeping only entries whose (run,entry) key is NOT in keys2_set
    auto df_diff = df_ok.Filter(
        [&keys2_set](const ActRoot::MergerData& m)
        {
            uint64_t key = (uint64_t(m.fRun) << 32) | uint64_t(uint32_t(m.fEntry));
            return keys2_set.find(key) == keys2_set.end();
        },
        {"MergerData"});

    auto count_diff = *df_diff.Count();
    std::cout << "Entries with fFlag == 'ok' in df but not in df2: " << count_diff << std::endl;

    // Get the events
    std::ofstream outFile("./Outputs/different_ok_events_preANDpost_Task0.dat");
    df_diff.Foreach(
        [&outFile](const ActRoot::MergerData& m)
        {
            m.Stream(outFile);
        },
        {"MergerData"});
    outFile.close();
}