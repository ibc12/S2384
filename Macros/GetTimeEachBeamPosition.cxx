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
void GetTimeEachBeamPosition()
{
    bool separateL1threshold = true;
    std::string which = "after";
    const int thresholdRun = 58;

    if(separateL1threshold && which != "before" && which != "after")
    {
        std::cout << "ERROR: 'which' debe ser \"before\" o \"after\" (recibido: " << which << ")\n";
        return;
    }

    // Selección de runs según flags
    auto isSelected = [&](int run)
    {
        if(!separateL1threshold)
            return true;
        if(run == thresholdRun)
            return false; // el run 58 nunca se incluye
        if(which == "before")
            return run < thresholdRun;
        return run > thresholdRun; // "after"
    };

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
    auto hRuns = gated.Histo1D({"hRuns", "Beam counts per run;Run;Counts", 122 - 19 + 1, 18.5, 122.5}, "fRun");

    // ---------------------------------------------------
    // 4. Sumar cuentas por periodo (solo runs seleccionados)
    // ---------------------------------------------------
    double total = 0;

    double pre7Li = 0;              // 19–65
    double post7Li_preMeasure = 0;  // 95–107
    double post7Li_postMeasure = 0; // 108–122

    for(int bin = 1; bin <= hRuns->GetNbinsX(); ++bin)
    {
        int run = static_cast<int>(hRuns->GetBinCenter(bin));
        double counts = hRuns->GetBinContent(bin);

        if(counts == 0 || !isSelected(run))
            continue;

        total += counts;

        if(run >= 19 && run <= 65)
            pre7Li += counts;
        else if(run >= 95 && run <= 107)
            post7Li_preMeasure += counts;
        else if(run >= 108 && run <= 122)
            post7Li_postMeasure += counts;
    }

    if(total == 0)
    {
        std::cout << "ERROR: no hay cuentas en los runs seleccionados\n";
        return;
    }

    // ---------------------------------------------------
    // 5. Escribir output .dat (nombre según flags)
    // ---------------------------------------------------
    std::string suffix = "";
    if(separateL1threshold)
        suffix = "_" + which + "L1threshold"; // _beforeL1threshold / _afterL1threshold

    std::string outname = "./Outputs/beamEmittancePeriods_11Li" + suffix + ".dat";
    std::ofstream out(outname);

    out << "pre7Li " << pre7Li / total << "\n";
    out << "post7Li_preMeasure " << post7Li_preMeasure / total << "\n";
    out << "post7Li_postMeasure " << post7Li_postMeasure / total << "\n";

    out.close();

    std::cout << "Output written to " << outname << "\n";

    // ---------------------------------------------------
    // 6. Plot de control
    // ---------------------------------------------------
    auto* c1 = new TCanvas("c1", "Beam counts per run");
    hRuns->DrawClone();
}