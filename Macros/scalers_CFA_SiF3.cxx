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

void scalers_CFA_SiF3()
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