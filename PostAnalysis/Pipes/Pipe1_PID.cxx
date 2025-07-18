#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TH2.h"
#include "TString.h"

#include <map>
#include <string>

void Pipe1_PID(std::string beam, std::string target, std::string light)
{
    // Read data
    ActRoot::DataManager dataman {"../configs/data.conf", ActRoot::ModeType::EMerge};
    auto chain {dataman.GetChain()};

    // RDataFrame
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {*chain};

    // LIGHT particle
    // Define lambda functions
    // 1-> Stopped in first silicon layer
    auto lambdaOne {[](ActRoot::MergerData& m) { return m.fLight.GetNLayers() == 1; }};
    // 2-> In two layers
    auto lambdaTwo {[](ActRoot::MergerData& m)
                    {
                        if(m.fLight.GetNLayers() == 2)
                            return (m.fLight.GetLayer(0) == "f0" && m.fLight.GetLayer(1) == "f1");
                        else
                            return false;
                    }};

    // Fill histograms
    std::map<std::string, ROOT::TThreadedObject<TH2D>> hsgas, hstwo;
    // Histogram models
    auto hGasSil {new TH2D {"hGasSil", ";E_{Sil} [MeV];#Delta E_{gas} [arb. units]", 300, 0, 30, 800, 0, 4000}};
    auto hTwoSils {new TH2D {"hTwoSils", ";#DeltaE_{0} [MeV];#DeltaE_{1} [MeV]", 300, 0, 30, 300, 0, 30}};
    for(const auto& layer : {"f0", "l0", "r0"})
    {
        hsgas.emplace(layer, *hGasSil);
        hsgas[layer]->SetTitle(TString::Format("%s", layer));
    }
    hstwo.emplace("f0-f1", *hTwoSils);
    hstwo["f0-f1"]->SetTitle("f0-f1");

    // Fill them
    df.Foreach(
        [&](ActRoot::MergerData& m)
        {
            if(lambdaOne(m)) // Gas-E0 PID
            {
                auto layer {m.fLight.GetLayer(0)};
                if(hsgas.count(layer))
                    hsgas[layer]->Fill(m.fLight.fEs.front(), m.fLight.fQave);
            }
            else if(lambdaTwo(m)) // E0-E1 PID
            {
                hstwo["f0-f1"]->Fill(m.fLight.fEs[0], m.fLight.fEs[1]);
            }
        },
        {"MergerData"});

    // If cuts are present, apply them
    ActRoot::CutsManager<std::string> cuts;
    // Gas PID
    cuts.ReadCut("gas", TString::Format("./Cuts/pid_%s_%s_gas.root", target.c_str(), light.c_str()).Data());
    // Two sils PID
    cuts.ReadCut("sils", TString::Format("./Cuts/pid_%s_%s_sils.root", target.c_str(), light.c_str()).Data());
    if(cuts.GetCut("gas") || cuts.GetCut("sils"))
    {
        // Apply PID and save in file
        auto gated {df.Filter(
            [&](ActRoot::MergerData& m)
            {
                if(cuts.GetCut("gas") && lambdaOne(m))
                    return cuts.IsInside("gas", m.fLight.fEs[0], m.fLight.fQave);
                else if(cuts.GetCut("sils") && lambdaTwo(m))
                    return cuts.IsInside("sils", m.fLight.fEs[0], m.fLight.fEs[1]);
                else
                    return false;
            },
            {"MergerData"})};
        auto name {TString::Format("./Outputs/tree_pid_%s_%s.root", target.c_str(), light.c_str())};
        std::cout << "Saving PID_Tree in file : " << name << '\n';
        gated.Snapshot("PID_Tree", name.Data());
    }

    // Draw
    auto* c0 {new TCanvas {"c0", "PID canvas"}};
    c0->DivideSquare(4);
    int p {1};
    for(auto& [layer, h] : hsgas)
    {
        c0->cd(p);
        h.Merge()->DrawClone("colz");
        p++;
    }
    for(auto& [layer, h] : hstwo)
    {
        c0->cd(p);
        h.Merge()->DrawClone("colz");
        p++;
    }
}
