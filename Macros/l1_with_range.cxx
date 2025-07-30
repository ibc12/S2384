#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActTypes.h"

#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TROOT.h"
void l1_with_range()
{
    ActRoot::DataManager dataman {"../configs/data.conf", ActRoot::ModeType::EMerge};
    dataman.SetRuns(55, 64);
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(chain2.get());

    // Dataframe
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {*chain};

    auto gated {
        df.Filter(
              [](ActRoot::MergerData& mer, ActRoot::ModularData& mod)
              {
                  if(mod.Get("GATCONF") != 8)
                      return false;
                  if(mer.fLightIdx == -1)
                      return false;
                  return true;
              },
              {"MergerData", "ModularData"})
            .Define("range", [](ActRoot::MergerData& mer) { return (mer.fBraggP - mer.fRP).R(); }, {"MergerData"})
            .Define("qtotal", [](ActRoot::MergerData& mer) { return mer.fLight.fQtotal; }, {"MergerData"})
            .Define("tl", [](ActRoot::MergerData& mer) { return mer.fLight.fTL; }, {"MergerData"})};

    // Book histograms
    auto hL1PID {
        gated.Histo2D({"hL1PID", "L1 PID;Range [mm];Qtotal [au]", 600, 0, 300, 2000, 0, 200000}, "range", "qtotal")};
    auto hL1TL {gated.Histo2D({"hL1TL", "L1 PID;TL [mm];Qtotal [au]", 600, 0, 300, 2000, 0, 200000}, "tl", "qtotal")};

    // Plot
    auto* c0 {new TCanvas {"c0", "l1 plot"}};
    c0->DivideSquare(4);
    c0->cd(1);
    hL1PID->DrawClone("colz");
    c0->cd(2);
    hL1TL->DrawClone("colz");
}
