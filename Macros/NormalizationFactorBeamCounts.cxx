#include "ActDataManager.h"
#include "ActSilData.h"
#include "ActTPCData.h"
#include "ActModularData.h"
#include "ActModularData.h"
#include "ActMergerData.h"
#include "ActTPCData.h"

#include <ROOT/RDataFrame.hxx>

#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"
#include "TLegend.h"

void NormalizationFactorBeamCounts()
{
    // === 1. Volvemos a cargar los datos ===
    ActRoot::DataManager dataman{"../configs/data_all.conf", ActRoot::ModeType::EReadSilMod};
    auto chain{dataman.GetChain()};

    auto friend1{dataman.GetChain(ActRoot::ModeType::EMerge)};
    chain->AddFriend(friend1.get());

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df{*chain};

    // === 2. Filtro CFA ===
    auto gatedCFA = df.Filter([](ActRoot::ModularData &d){ return d.Get("GATCONF") == 64; }, {"ModularData"});

    // Rango de runs
    auto minRun = *df.Min<int>("fRun");
    auto maxRun = *df.Max<int>("fRun");
    int nbinsX = maxRun - minRun + 1;

    // Histograma CFA counts por run
    auto cfa_counts = gatedCFA.Histo1D({"cfaCounts", "CFA counts per run;Run;Counts",
                                        nbinsX, minRun - 0.5, maxRun + 0.5}, "fRun");

    // Clonamos histograma para modificarlo
    auto cfa_h = (TH1D*)cfa_counts->Clone("cfa_h_scaled");
    cfa_h->Scale(300.0); // === 3. Multiplicamos por el factor 300 ===

    // === 4. Dividir bin a bin por el valor del TGraph guardado ===
    // Cargamos el TGraph desde archivo (ajusta el path/nombre según lo guardaste)
    TFile fIn("./Outputs/CFA_F2_ratio_normalized_by_first_run.root"); // <-- tu archivo con el TGraph
    auto gRatio = (TGraph*)fIn.Get("Graph"); // <-- nombre del TGraph
    if (!gRatio) {
        std::cerr << "Error: no se pudo cargar gRatio del archivo ratioGraph.root" << std::endl;
        return;
    }

    double totalCounts = 0; // acumulador

    for (int i = 1; i <= cfa_h->GetNbinsX(); i++) {
        double run = cfa_h->GetBinCenter(i);
        double ratioVal = gRatio->Eval(run); // interpola en el run correspondiente
        if (ratioVal != 0) {
            double newVal = cfa_h->GetBinContent(i) / ratioVal;
            cfa_h->SetBinContent(i, newVal);
            totalCounts += newVal; // acumulamos
        } else {
            cfa_h->SetBinContent(i, 0); // evitar div/0
        }
    }

    std::cout << "========================================" << std::endl;
    std::cout << "   Total CFA*300 / (CFA/F2 ratio) = " << totalCounts << std::endl;
    std::cout << "========================================" << std::endl;

    // === 5. Dibujar ===
    auto c1 = new TCanvas("c1", "CFA scaled / Ratio", 1000, 600);
    cfa_h->SetTitle("CFA*300 / (CFA/F2 ratio);Run;Scaled counts");
    cfa_h->SetMarkerStyle(20);
    cfa_h->SetMarkerColor(kRed+1);
    cfa_h->Draw("hist");
}