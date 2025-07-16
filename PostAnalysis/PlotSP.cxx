#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActTypes.h"

#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TH2.h"
#include "TROOT.h"
#include "TString.h"

#include <map>
#include <string>

void PlotSP()
{
    // Read the data using the data.conf file
    ActRoot::DataManager dataman {"../configs/data.conf", ActRoot::ModeType::EMerge};
    auto chain {dataman.GetChain()}; // Get all Merge files for Runs in a single TChain
    // Add friends if necessary
    auto friend1 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(friend1.get());

    // Build the RDataFrame
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {*chain};

    // Gate on events
    auto gated {df.Filter(
        [](ActRoot::MergerData& m)
        {
            if(!m.fLight.HasSP())
                return false;
            return true;
        },
        {"MergerData"})};

    // Fill histograms
    int nsils {11};
    std::map<std::string, std::map<int, ROOT::TThreadedObject<TH2D>>> hs;
    // Histogram model
    auto* h2d {new TH2D {"h2d", "SP;X or Y [pad];Z [btb]", 300, 0, 300, 300, 0, 300}};
    for(const auto& layer : {"f0", "l0", "r0"})
    {
        for(int s = 0; s < nsils; s++)
        {
            hs[layer].emplace(s, *h2d);
            hs[layer][s]->SetTitle(TString::Format("Layer %s", layer));
        }
    }
    gated.Foreach(
        [&](ActRoot::MergerData& m)
        {
            // No need to check for SP, it has already been done in Filter
            auto layer {m.fLight.fLayers.front()}; // ensured to have at least size >= 1
            auto n {m.fLight.fNs.front()};
            auto sp {m.fLight.fSP};
            if(hs.count(layer))
            {
                if(hs[layer].count(n))
                {
                    if(layer == "f0")
                        hs[layer][n].Get()->Fill(sp.Y(), sp.Z());
                    else
                        hs[layer][n].Get()->Fill(sp.X(), sp.Z());
                }
            }
        },
        {"MergerData"});


    // Draw
    auto* c0 {new TCanvas {"c0", "SP canvas"}};
    c0->DivideSquare(4);
    int p {1};
    for(auto& [layer, hsils] : hs)
    {
        c0->cd(p);
        int idx {};
        for(auto& [s, h] : hsils)
        {
            auto color {idx + 1};
            if(color == 10) // 10 is white, as well as 0
                color = 46;
            auto opts {(idx == 0) ? "scat" : "scat same"};
            // Merge histos from threads
            h.Merge();
            // Set color
            h.GetAtSlot(0)->SetMarkerColor(color);
            // Draw
            h.GetAtSlot(0)->Draw(opts);
            idx++;
        }
        p++;
    }
}
