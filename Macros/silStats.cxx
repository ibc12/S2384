#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TH2D.h"
#include "TString.h"

void silStats()
{
    ActRoot::DataManager dataman {"../configs/data.conf", ActRoot::ModeType::EMerge};
    auto chain {dataman.GetChain()};
    ROOT::RDataFrame df {*chain};

    // Build the map
    std::map<std::string, TH2D*> hs;
    for(const auto& layer : {"l0", "r0", "f0", "f2"})
    {
        hs[layer] = new TH2D(TString::Format("h%s", layer), TString::Format("Sil stats for %s;Col;Row", layer), 4, 0, 4,
                             4, 0, 4);
    }

    // Fill!
    df.Foreach(
        [&](ActRoot::MergerData& mer)
        {
            if(mer.fLight.GetNLayers() == 1)
            {
                auto layer {mer.fLight.GetLayer(0)};
                if(!hs.count(layer))
                    return;
                auto n {mer.fLight.fNs[0]};
                auto e {mer.fLight.fEs[0]};
                int row {};
                int col {};
                if(layer == "f0")
                {
                    row = n / 4;
                    col = n % 4;
                }
                else if(layer == "l0" || layer == "r0")
                {
                    row = n / 3;
                    col = n % 3;
                }
                else
                    return;
                hs[layer]->Fill(col, row);
            }
            if(mer.fHeavy.GetNLayers() > 0)
            {
                if(mer.fHeavy.fLayers[0] == "f2")
                {
                    auto n {mer.fHeavy.fNs[0]};
                    hs["f2"]->Fill(static_cast<int>(n % 2), static_cast<int>(n / 2));
                }
            }
        },
        {"MergerData"});

    // Draw
    auto* c0 {new TCanvas {"c0", "Sil stats canvas"}};
    c0->DivideSquare(hs.size());
    int p {1};
    for(auto& [layer, h] : hs)
    {
        c0->cd(p);
        h->Draw("colz");
        p++;
    }
}
