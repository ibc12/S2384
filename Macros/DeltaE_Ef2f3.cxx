#include "ActDataManager.h"
#include "ActSilData.h"
#include "ActTPCData.h"
#include "ActModularData.h"

#include <ROOT/RDataFrame.hxx>

#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"

void DeltaE_Ef2f3()
{
    ActRoot::DataManager dataman {"../configs/data.conf", ActRoot::ModeType::EReadSilMod};
    auto chain {dataman.GetChain()};

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {*chain};

    df.Describe().Print();

    // Define map of histograms
    std::map<int, TH2D*> hs;
    for(int s = 0; s < 4; s++)
    {
        hs[s] = new TH2D(Form("hDeltaEf2f3_s%d", s), Form("#DeltaE-E f2-f3 for %d;E f3 [MeV];E f2 [MeV]", s), 450, 0,
                         80, 160, 0, 25);
    }

    auto dfFilterGATCONF = df.Filter(
        [](ActRoot::ModularData& d) { return d.Get("GATCONF") == 4; }, {"ModularData"});

    auto gated {
        dfFilterGATCONF.Filter([](ActRoot::SilData& d) { return (d.fSiE["f2"].size() && d.fSiE["f3"].size()); }, {"SilData"})};

    gated.Foreach(
        [&](ActRoot::SilData& d)
        {
            for(int i = 0; i < d.fSiE["f2"].size(); i++)
            {
                auto n {d.fSiN["f2"][i]};
                auto e2 {d.fSiE["f2"][i]};
                auto e3 {d.fSiE["f3"].front()};
                if(e2 < 0.5)
                    continue;
                hs[n]->Fill(e3, e2);
            }
        },
        {"SilData"});
    //auto dfnew = df.Define("deltaEf2",
    //                       [](ActRoot::SilData& d)
    //                       {
    //                           int pad {0};
    //                           if(!(d.fSiN["f2"].empty()) && !(d.fSiN["f3"].empty()))
    //                           {
    //                               auto& f2 {d.fSiN["f2"]};
    //                               auto it {std::find(f2.begin(), f2.end(), pad)};
    //                               if(it != f2.end())
    //                               {
    //                                   auto idx {std::distance(f2.begin(), it)};
    //                                   auto e {d.fSiE["f2"][idx]};
    //                                   if(e > 0.5f)
    //                                       return e;
    //                                   else
    //                                       return -1.0f;
    //                               }
    //                               else
    //                                   return -1.0f;
    //                           }
    //                           else
    //                               return -1.0f;
    //                       },
    //                       {"SilData"})
    //                 .Define("Ef3",
    //                         [](ActRoot::SilData& d)
    //                         {
    //                             int pad {0};
    //                             if(!(d.fSiN["f2"].empty()) && !(d.fSiN["f3"].empty()))
    //                             {
    //                                 auto& f3 {d.fSiN["f3"]};
    //                                 auto it {std::find(f3.begin(), f3.end(), pad)};
    //                                 if(it != f3.end())
    //                                 {
    //                                     auto idx {std::distance(f3.begin(), it)};
    //                                     auto e {d.fSiE["f3"][idx]};
    //                                     if(e > 0.5f)
    //                                         return e;
    //                                     else
    //                                         return -1.0f;
    //                                 }
    //                                 else
    //                                     return -1.0f;
    //                             }
    //                             else
    //                                 return -1.0f;
    //                         },
    //                         {"SilData"});
    // auto histo {dfnew.Histo2D({"hDeltaEf2f3", "Delta E f2 vs f3;E f3 [MeV];E f2 [MeV]", 450, 0, 80, 160, 0, 25}, "Ef3",                           "deltaEf2")};

    TCanvas* c1 = new TCanvas("c1", "Delta E f2 vs f3", 800, 600);
    c1->DivideSquare(hs.size());
    for(int s = 0; s < hs.size(); s++)
    {
        c1->cd(s + 1);
        hs[s]->SetStats(0);
        hs[s]->DrawClone("colz");
    }

    auto* c2 = new TCanvas("c2", "Delta E f2 vs f3 for 1 cuadrant sils", 800, 600);
    c2->cd();
    hs[1]->DrawClone("colz");
       
}
