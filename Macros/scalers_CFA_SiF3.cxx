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

void do_it_respect_to_time()
{
    // Get Sil data
    ActRoot::DataManager dataman{"../configs/data.conf", ActRoot::ModeType::EReadSilMod};
    auto chain{dataman.GetChain()};

    // Build the RDataFrame
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df{*chain};

    // Now print the results
    // df.Describe().Print();

    // Define CFA/SiF3 scalers division
    auto df2 = df.Define("CFA_SiF3", [](ActRoot::ModularData &d)
                         {
                             double cfa{d.Get("SCA_CFA")};
                             double sif3{d.Get("SCA_SIF3")};
                             if (sif3 != 0)
                                 return cfa / sif3;
                             else
                                 return -1.0; // or some other value indicating invalid division
                         },
                         {"ModularData"})
                   .Define("time", [](ActRoot::ModularData &d)
                           {
                            auto timeh_up {d.Get("CTR_TIMEH_UP")};
                            auto timeh {d.Get("CTR_TIMEH")};
                            auto timeml_up {d.Get("CTR_TIMEML_UP")};
                            auto timeml {d.Get("CTR_TIMEML")};

                            return (timeml + timeml_up * 2e16 + timeh * 2e32 + timeh_up * 2e48) * 12.5e-9; }, {"ModularData"}); // suposing 12.5 ns per second

    auto max_time = *df2.Max<double>("time");
    std::cout << "Tiempo máximo: " << max_time << std::endl;

    // Plot the value
    auto histo = df2.Histo1D({"hCFA_SiF3", "CFA/SiF3 Ratio;CFA/SiF3;Counts", 100, -2, 20}, "CFA_SiF3");
    auto histo2 = df2.Histo2D({"hCFA_SiF3_vs_time", "CFA/SiF3 vs Time;Time [arb. units];CFA/SiF3", 10000, 0, 5.25e26, 100, -2, 20}, "time", "CFA_SiF3");
    // Create canvas and draw histogram
    TCanvas *c1 = new TCanvas("c1", "CFA/SiF3 Ratio", 800, 600);
    c1->Divide(1,2);
    c1->cd(1);
    histo->DrawClone();
    c1->cd(2);
    histo2->DrawClone("colz");
}

void scalers_CFA_SiF3()
{
    // Get Sil data
    ActRoot::DataManager dataman{"../configs/data.conf", ActRoot::ModeType::EReadSilMod};
    auto chain{dataman.GetChain()};

    // Build the RDataFrame
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df{*chain};

    // Define CFA/SiF3 ratio and entry index
    auto df2 = df.Define("CFA_SiF3", [](ActRoot::ModularData &d) {
                            double cfa{d.Get("SCA_CFA")};
                            double sif3{d.Get("SCA_SIF3")};
                            if (sif3 != 0)
                                return cfa / sif3;
                            else
                                return -1.0; // marcador de división inválida
                        },
                        {"ModularData"})
                   .Define("entry", [](ULong64_t i) {
                            return static_cast<double>(i);
                        }, {"rdfentry_"}); // índice del evento en el dataframe

    // Número total de eventos
    auto nentries = *df2.Count();
    std::cout << "Número total de eventos: " << nentries << std::endl;

    // Elegir número de bins en X (máximo 5000 para no saturar)
    int maxBinsX = 5000;
    int nbinsX = (nentries < static_cast<ULong64_t>(maxBinsX)) ? static_cast<int>(nentries) : maxBinsX;

    // Histogramas
    auto histo = df2.Histo1D({"hCFA_SiF3",
                              "CFA/SiF3 Ratio;CFA/SiF3;Counts",
                              100, -2, 20},
                             "CFA_SiF3");

    auto histo2 = df2.Histo2D({"hCFA_SiF3_vs_entry",
                               "CFA/SiF3 vs Entry;Entry index;CFA/SiF3",
                               nbinsX, 0, static_cast<double>(nentries),
                               100, -2, 20},
                              "entry", "CFA_SiF3");

    // Clonar el TH2D para aplicar el filtro
    auto hFiltered = (TH2D*)histo2->Clone("hFiltered");

    // Recorrer bins y dejar solo los >60
    for (int ix = 1; ix <= hFiltered->GetNbinsX(); ix++) {
        for (int iy = 1; iy <= hFiltered->GetNbinsY(); iy++) {
            if (hFiltered->GetBinContent(ix, iy) <= 60) {
                hFiltered->SetBinContent(ix, iy, 0);
            }
        }
    }

    // Canvas
    auto c1 = new TCanvas("c1", "CFA/SiF3 Analysis", 1200, 800);
    c1->Divide(1,2);
    c1->cd(1);
    histo->DrawClone();
    c1->cd(2);
    hFiltered->Draw("colz");
}
