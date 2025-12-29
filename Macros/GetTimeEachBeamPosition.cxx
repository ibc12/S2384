#include "ActDataManager.h"
#include "ActModularData.h"
#include "ActTPCData.h"

#include <ROOT/RDataFrame.hxx>

#include <TCanvas.h>
#include <TH1D.h>

#include <fstream>
#include <iostream>

void GetTimeEachBeamPosition()
{
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
    // 4. Sumar cuentas por periodo
    // ---------------------------------------------------
    double total = 0;

    double pre7Li = 0;              // 19–65
    double post7Li_preMeasure = 0;  // 95–107
    double post7Li_postMeasure = 0; // 108–122

    for(int bin = 1; bin <= hRuns->GetNbinsX(); ++bin)
    {
        int run = static_cast<int>(hRuns->GetBinCenter(bin));
        double counts = hRuns->GetBinContent(bin);

        if(counts == 0)
            continue;

        total += counts;

        if(run >= 19 && run <= 65)
            pre7Li += counts;
        else if(run >= 95 && run <= 107)
            post7Li_preMeasure += counts;
        else if(run >= 108 && run <= 122)
            post7Li_postMeasure += counts;
    }

    // ---------------------------------------------------
    // 5. Escribir output .dat (formato mínimo)
    // ---------------------------------------------------
    std::string outname = "./Outputs/beamEmittancePeriods_11Li.dat";
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
