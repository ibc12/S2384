#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>

#include <map>
#include <set>
#include <fstream>
#include <string>

#include "ActDataManager.h"
#include "ActModularData.h"
#include "ActTPCData.h"
#include "ActMergerData.h"

void GetWorkingTimeForEachSilicon() {
    // Get Sil data
    ActRoot::DataManager dataman{"../configs/data.conf", ActRoot::ModeType::EReadSilMod};
    auto chain{dataman.GetChain()};
    auto friend1{dataman.GetChain(ActRoot::ModeType::EFilter)};
    chain->AddFriend(friend1.get());
    auto friend2{dataman.GetChain(ActRoot::ModeType::EMerge)};
    chain->AddFriend(friend2.get());

    // Build the RDataFrame
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df{*chain};

    // ---------------------------------------------------
    // 1. Filtrar a eventos de haz (tu "gated")
    // ---------------------------------------------------
    auto gated = df.Filter(
                      [](ActRoot::ModularData &d) { return d.Get("GATCONF") == 64; },
                      {"ModularData"})
                      .Filter(
                          [](ActRoot::TPCData &d) { return d.fClusters.size() >= 1; },
                          {"TPCData"})
                      .Filter(
                          [](ActRoot::TPCData &d) {
                              auto [xmin, xmax] = d.fClusters.front().GetXRange();
                              return (xmin < 5 && xmax > 120);
                          },
                          {"TPCData"});
    gated.Describe().Print();
    // ---------------------------------------------------
    // 2. Histograma de cuentas de haz por run
    // ---------------------------------------------------
    auto hRuns = gated.Histo1D(
        {"hRuns", "Beam counts per run", 122-19+1, 18.5, 122.5}, "fRun");

    std::map<int, double> beamCounts;
    for (int i = 1; i <= hRuns->GetNbinsX(); i++) {
        int run = i; // bin number = run (si tu numeración de run es así)
        double val = hRuns->GetBinContent(i);
        if (val > 0) {
            beamCounts[run] = val;
        }
    }

    // ---------------------------------------------------
    // 3. Runs malos por detector
    // ---------------------------------------------------
    std::map<std::string, std::set<int>> badRuns;
    badRuns["l0_9"] = {}; // siempre malo → lo tratamos después
    badRuns["f0_2"] = {};
    badRuns["r0_3"] = {};

    // Casos especiales
    badRuns["r0_3"].insert({30, 31, 33});
    badRuns["r0_5"].insert({30, 31, 33});
    badRuns["r0_all"].insert({32});
    badRuns["r0_all"].insert({34});
    badRuns["l0_all"].insert({34});
    for (int r = 36; r < 45; r++) {
        badRuns["f0_5"].insert(r);
    }
    badRuns["r0_2"].insert({116, 117});

    // ---------------------------------------------------
    // 4. Calcular cuentas buenas por detector
    // ---------------------------------------------------
    std::ofstream out("detector_counts.txt");

    // Total de cuentas en todos los runs
    double totalAll = 0;
    for (auto &[run, c] : beamCounts) {
        totalAll += c;
    }

    for (auto det : {"f0", "r0", "l0"}) {
        for (int i = 0; i <= 11; i++) {
            std::string name = Form("%s_%d", det, i);

            double totalGood = 0;

            for (auto &[run, c] : beamCounts) {
                bool isBad = false;

                // Detector siempre malo
                if ((det == std::string("l0") && i == 9) ||
                    (det == std::string("f0") && i == 2) ||
                    (det == std::string("r0") && i == 3)) {
                    isBad = true;
                }

                // Runs listados como malos
                if (badRuns[name].count(run)) {
                    isBad = true;
                }

                // Runs donde todo el lado está apagado
                if (badRuns[std::string(det) + "_all"].count(run)) {
                    isBad = true;
                }

                if (!isBad) {
                    totalGood += c;
                }
            }

            // Normalización en base 1
            double goodFrac = (totalAll > 0) ? (totalGood / totalAll) : 0.0;
            out << name << ": " << goodFrac << "\n";
        }
    }

    out.close();
    

    // Plot histogram of counts per run
    auto *c1 = new TCanvas("c1", "Beam counts per run");
    hRuns->DrawClone();
}
