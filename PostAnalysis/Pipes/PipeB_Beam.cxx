#include "ActDataManager.h"
#include "ActModularData.h"
#include "ActMergerData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TROOT.h"
#include "TString.h"

#include <atomic>
#include <fstream>
#include <stdexcept>

void PipeB_Beam(const std::string& beam)
{
    std::string dataconf {};
    if(beam == "11Li")
        dataconf = "./../configs/data_11Li.conf";
    else if(beam == "7Li")
        dataconf = "./../configs/data_7Li.conf";
    else
        throw std::runtime_error("Beam cannot differ from 11Li or 7Li");

    ROOT::EnableImplicitMT();
    // Read data
    ActRoot::DataManager datman {dataconf, ActRoot::ModeType::EReadSilMod};
    auto chain {datman.GetChain()};
    auto chain2 {datman.GetChain(ActRoot::ModeType::EMerge)};
    chain->AddFriend(chain2.get());
    ROOT::RDataFrame df {*chain};

    // Get GATCONF values y el Run
    auto def {df.Define("GATCONF", [](ActRoot::ModularData& mod) { return static_cast<int>(mod.fLeaves["GATCONF"]); },
                        {"ModularData"})
                  .Define("Run", [](ActRoot::MergerData& m) { return m.fRun; }, {"MergerData"})};

    // Book histograms
    auto hGATCONF {def.Histo1D("GATCONF")};

    // CFA counters per run, one per slot to avoid missing counts due to parallel execution
    auto nSlots {ROOT::IsImplicitMTEnabled() ? ROOT::GetThreadPoolSize() : 1u};
    if(nSlots == 0)
        nSlots = 1;
    std::vector<std::map<int, unsigned long>> cfaPerRunSlots(nSlots);

    def.ForeachSlot(
        [&](unsigned int slot, int gatconf, int run)
        {
            if(gatconf == 64)
                cfaPerRunSlots[slot][run]++;
        },
        {"GATCONF", "Run"});

    // Fuse the counts from all slots into a single map
    std::map<int, unsigned long> cfaPerRun {};
    for(const auto& slotMap : cfaPerRunSlots)
        for(const auto& [run, counts] : slotMap)
            cfaPerRun[run] += counts;

    // Total, for the report
    unsigned long cfa {};
    for(const auto& [run, counts] : cfaPerRun)
        cfa += counts;

    // Draw
    auto* c0 {new TCanvas {"c00", "Pipe 0 canvas 0"}};
    hGATCONF->DrawClone();

    // Print report
    std::cout << "===== GATCONF report =====" << '\n';
    std::cout << "-> CFA/div = " << cfa << '\n';
    std::cout << "==========================" << '\n';

    // save CFA counts per run, ordered from highest to lowest run
    auto filename {TString::Format("../Fits/norm/cfa_perRun_%s.dat", beam.c_str())};
    std::ofstream streamer {filename};
    if(!streamer)
        throw std::runtime_error("No se pudo abrir el fichero " + filename);
    for(auto it = cfaPerRun.rbegin(); it != cfaPerRun.rend(); ++it)
        streamer << it->first << "  " << it->second << '\n';
    streamer.close();

    std::cout << "-> Fichero escrito: " << filename << '\n';
}