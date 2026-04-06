#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"

#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <mutex>
#include <string>
#include <unordered_set>
#include <vector>

void DiffrentCountsWithTask0()
{
    std::string treeName = "ACTAR_Merged";

    // 👉 Múltiples runs (ajusta rango)
    std::vector<std::string> files;
    std::vector<std::string> files2;

    std::vector<int> runs;

    for(int r = 19; r <= 122; ++r)
    {
        // excluir runs
        if(r == 25 || r == 35 || r == 58 || r == 68 || (r >= 83 && r <= 94))
            continue;

        runs.push_back(r);
    }

    for(int run : runs)
    {
        files.push_back(Form("../../RootFiles/Merger/Merged_Run_%04d.root", run));
        files2.push_back(
            Form("../../RootFiles/Merger_all_postMultiActionChange_RebinZ_preTask0/Merged_Run_%04d.root", run));
    }

    ROOT::RDataFrame df(treeName, files);
    ROOT::RDataFrame df2(treeName, files2);

    auto total = *df.Count();
    auto total2 = *df2.Count();

    std::cout << "Total entries df: " << total << std::endl;
    std::cout << "Total entries df2: " << total2 << std::endl;

    // ✅ Filtro correcto
    auto df_ok = df.Filter([](const ActRoot::MergerData& m) { return m.fFlag == "ok"; }, {"MergerData"});

    auto df2_ok = df2.Filter([](const ActRoot::MergerData& m) { return m.fFlag == "ok"; }, {"MergerData"});

    auto count_ok = *df_ok.Count();
    auto count2_ok = *df2_ok.Count();

    std::cout << "Entries OK df: " << count_ok << std::endl;
    std::cout << "Entries OK df2: " << count2_ok << std::endl;

    // ============================================================
    // 🔹 Construir sets de claves (run, entry)
    // ============================================================

    std::vector<uint64_t> keys1, keys2;
    std::mutex mtx1, mtx2;

    df_ok.Foreach(
        [&keys1, &mtx1](const ActRoot::MergerData& m)
        {
            uint64_t key = (uint64_t(m.fRun) << 32) | uint64_t(uint32_t(m.fEntry));
            std::lock_guard<std::mutex> lock(mtx1);
            keys1.push_back(key);
        },
        {"MergerData"});

    df2_ok.Foreach(
        [&keys2, &mtx2](const ActRoot::MergerData& m)
        {
            uint64_t key = (uint64_t(m.fRun) << 32) | uint64_t(uint32_t(m.fEntry));
            std::lock_guard<std::mutex> lock(mtx2);
            keys2.push_back(key);
        },
        {"MergerData"});

    std::unordered_set<uint64_t> set1(keys1.begin(), keys1.end());
    std::unordered_set<uint64_t> set2(keys2.begin(), keys2.end());

    std::cout << "Unique df_ok: " << set1.size() << std::endl;
    std::cout << "Unique df2_ok: " << set2.size() << std::endl;

    // ============================================================
    // 🔹 Diferencias en ambas direcciones
    // ============================================================

    auto df_only = df_ok.Filter(
        [&set2](const ActRoot::MergerData& m)
        {
            uint64_t key = (uint64_t(m.fRun) << 32) | uint64_t(uint32_t(m.fEntry));
            return set2.find(key) == set2.end();
        },
        {"MergerData"});

    auto df2_only = df2_ok.Filter(
        [&set1](const ActRoot::MergerData& m)
        {
            uint64_t key = (uint64_t(m.fRun) << 32) | uint64_t(uint32_t(m.fEntry));
            return set1.find(key) == set1.end();
        },
        {"MergerData"});

    auto count_df_only = *df_only.Count();
    auto count_df2_only = *df2_only.Count();

    std::cout << "In df but NOT in df2: " << count_df_only << std::endl;
    std::cout << "In df2 but NOT in df: " << count_df2_only << std::endl;

    // ============================================================
    // 🔹 Guardar resultados
    // ============================================================

    std::ofstream out1("./Outputs/different_ok_events_preANDpost_Task0.dat");

    df_only.Foreach([&out1](const ActRoot::MergerData& m) { m.Stream(out1); }, {"MergerData"});

    df2_only.Foreach([&out1](const ActRoot::MergerData& m) { m.Stream(out1); }, {"MergerData"});

    out1.close();
}