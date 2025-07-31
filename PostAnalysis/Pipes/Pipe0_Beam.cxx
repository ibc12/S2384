#include "ActDataManager.h"
#include "ActModularData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TROOT.h"

#include <atomic>
#include <stdexcept>

void Pipe0_Beam(const std::string& beam)
{
    std::string dataconf {};
    if(beam == "11Li")
        dataconf = "./../configs/data.conf";
    else if(beam == "7Li")
        dataconf = "./../configs/data_7Li.conf";
    else
        throw std::runtime_error("Beam cannot differ from 11Li or 7Li");

    ROOT::EnableImplicitMT();
    // Read data
    ActRoot::DataManager datman {dataconf, ActRoot::ModeType::EReadSilMod};
    auto chain {datman.GetJoinedData()};
    ROOT::RDataFrame df {*chain};

    // Get GATCONF values
    auto def {df.Define("GATCONF", [](ActRoot::ModularData& mod) { return static_cast<int>(mod.fLeaves["GATCONF"]); },
                        {"ModularData"})};

    // Book histograms
    auto hGATCONF {def.Histo1D("GATCONF")};

    // And cound CFA triggers
    std::atomic<unsigned long int> cfa {};
    def.Foreach(
        [&](const int& gatconf)
        {
            if(gatconf == 64)
                cfa++;
        },
        {"GATCONF"});

    // Draw
    auto* c0 {new TCanvas {"c00", "Pipe 0 canvas 0"}};
    hGATCONF->DrawClone();

    // Print report
    std::cout << "===== GATCONF report =====" << '\n';
    std::cout << "-> CFA/div = " << cfa << '\n';
    std::cout << "==========================" << '\n';
}
