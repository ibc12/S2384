#include "ActDataManager.h"
#include "ActModularData.h"
#include "ActTPCData.h"

#include <ROOT/RDataFrame.hxx>

#include <TCanvas.h>
#include <TH1D.h>

#include <fstream>
#include <iostream>
#include <string>

// separateL1threshold = false -> todos los runs (comportamiento original)
// separateL1threshold = true  -> se separa por el cambio de threshold del trigger (run 58, excluido)
//     which = "before" -> runs < 58
//     which = "after"  -> runs > 58
void GetTimeEachL1TriggerFor1LiPreChange()
{
    const int thresholdRun = 58;
    const int run11LiEnd = 65;

    // ---------------------------------------------------
    // 1. Leer datos de 11Li
    // ---------------------------------------------------
    std::string dataconf = "./../configs/data_11Li.conf";
    ActRoot::DataManager dataman {dataconf, ActRoot::ModeType::EReadSilMod};

    auto chain = dataman.GetChain();
    auto friend1 = dataman.GetChain(ActRoot::ModeType::EFilter);
    auto friend2 = dataman.GetChain(ActRoot::ModeType::EMerge);

    chain->AddFriend(friend1.get());
    chain->AddFriend(friend2.get());

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {*chain};

    // ---------------------------------------------------
    // 2. Gate simple de eventos de haz
    // ---------------------------------------------------
    auto gated = df.Filter([](ActRoot::ModularData& d) { return d.Get("GATCONF") == 64; }, {"ModularData"})
                     .Filter([](ActRoot::TPCData& d) { return d.fClusters.size() >= 1; }, {"TPCData"});

    // ---------------------------------------------------
    // 3. Histograma de cuentas por run
    // ---------------------------------------------------
    auto hRuns = gated.Histo1D({"hRuns", "Beam counts per run;Run;Counts", 65 - 19 + 1, 18.5, 65.5}, "fRun");

    // ---------------------------------------------------
    // 4. Sumar cuentas por periodo (solo runs seleccionados)
    // ---------------------------------------------------
    double total = 0;

    double preTriggerChange = 0;              // 19–57
    double postTriggerChange_pre7Li = 0;      // 59–65

    for(int bin = 1; bin <= hRuns->GetNbinsX(); ++bin)
    {
        int run = static_cast<int>(hRuns->GetBinCenter(bin));
        double counts = hRuns->GetBinContent(bin);

        total += counts;

        if(run >= 19 && run < thresholdRun)
            preTriggerChange += counts;
        else if(run >= thresholdRun + 1 && run <= run11LiEnd)
            postTriggerChange_pre7Li += counts;
    }

    if(total == 0)
    {
        std::cout << "ERROR: no hay cuentas en los runs seleccionados\n";
        return;
    }

    // ---------------------------------------------------
    // 5. Escribir output .dat
    // ---------------------------------------------------

    std::string outname = "./Outputs/triggerL1Periods_11Li.dat";
    std::ofstream out(outname);

    out << "preTriggerChange " << preTriggerChange / total << "\n";
    out << "postTriggerChange_pre7Li " << postTriggerChange_pre7Li / total << "\n";

    out.close();

    std::cout << "Output written to " << outname << "\n";

    // ---------------------------------------------------
    // 6. Plot de control
    // ---------------------------------------------------
    auto* c1 = new TCanvas("c1", "Beam counts per run");
    hRuns->DrawClone();
}