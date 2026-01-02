#include "ActDataManager.h"
#include "ActSilData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TString.h"
#include "TVirtualPad.h"

#include <map>
#include <string>
#include <vector>

void rawStats()
{
    ActRoot::DataManager dataman {"../configs/data.conf", ActRoot::ModeType::EReadSilMod};
    auto chain {dataman.GetChain()};

    ROOT::RDataFrame df {*chain};

    std::map<std::string, std::vector<TH1D*>> hs;
    int nsil {12};
    for(const auto& layer : {"f0", "l0", "r0"})
    {
        hs[layer] = {};
        for(int s = 0; s < nsil; s++)
        {
            hs[layer].push_back(new TH1D {TString::Format("h%s%d", layer, s),
                                          TString::Format("Energies %s_%d;E [MeV];Counts", layer, s), 300, 0, 80});
        }
    }
    // Fill!
    df.Foreach(
        [&](ActRoot::SilData& sil)
        {
            for(const auto& layer : sil.GetLayers())
            {
                if(!hs.count(layer))
                    continue;
                auto mult {sil.fSiE[layer].size()};
                for(int m = 0; m < mult; m++)
                {
                    auto n {sil.fSiN[layer][m]};
                    hs[layer][n]->Fill(sil.fSiE[layer][m]);
                }
            }
        },
        {"SilData"});

    // Draw
    for(auto& [layer, vh] : hs)
    {
        auto* c {new TCanvas {TString::Format("c%s", layer.c_str()), TString::Format("Raw sil %s", layer.c_str())}};
        c->DivideSquare(vh.size());
        int p {1};
        for(auto& h : vh)
        {
            c->cd(p);
            gPad->SetLogy();
            h->Draw();
            p++;
        }
    }
}
