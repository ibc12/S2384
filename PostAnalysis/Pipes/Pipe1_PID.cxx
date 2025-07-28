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
    ActRoot::DataManager dataman{"../configs/data.conf", ActRoot::ModeType::EMerge};
    auto chain{dataman.GetChain()};

    // RDataFrame
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame dforigin{*chain};

    // Filter some of the silicons
    auto dfFilterLeftSilicon = dforigin.Filter([](ActRoot::MergerData &m)
                              {
            if(!m.fLight.fLayers.empty() && m.fLight.fLayers.front() == "l0")
            {
                if(!m.fLight.fNs.empty() && m.fLight.fNs.front() == 9)
                {
                    return false;
                }
                else
                    return true;
            }
            else    
                return true; }, {"MergerData"});

    // Filter the right wall for the runs in 27/07 night
    auto df = dfFilterLeftSilicon.Filter([](ActRoot::MergerData &m) {
        if(m.fRun > 29 && m.fRun < 37)
        {
            if(m.fRun == 30 || m.fRun == 31)
            {
                if(!m.fLight.fLayers.empty() && m.fLight.fLayers.front() == "r0")
                {
                    if(!m.fLight.fNs.empty() && (m.fLight.fNs.front() == 3 || m.fLight.fNs.front() == 5))
                    {
                        return false;
                    }
                    else
                        return true;
                }
                else    
                    return true;
            }
            if(m.fRun == 32)
            {
                if(!m.fLight.fLayers.empty() && m.fLight.fLayers.front() == "r0")
                {
                    return false;
                }
                else
                    return true;
            }
            if(m.fRun == 34)
            {
                if(!m.fLight.fLayers.empty() && (m.fLight.fLayers.front() == "r0" || m.fLight.fLayers.front() == "l0"))
                {
                    return false;
                }
                else
                    return true;
            }
            else
                return true;

        }
        else
            return true;
    },
                   {"MergerData"});

    // print example of df
    // df.Describe().Print();

    // LIGHT particle
    // Define lambda functions
    // 1-> Stopped in first silicon layer
    auto lambdaOne{[](ActRoot::MergerData &m)
                   { return m.fLight.GetNLayers() == 1; }};
    // 2-> In two layers
    auto lambdaTwo{[](ActRoot::MergerData &m)
                   {
                       if (m.fLight.GetNLayers() == 2)
                           return (m.fLight.GetLayer(0) == "f0" && m.fLight.GetLayer(1) == "f1");
                       else
                           return false;
                   }};

    // Fill histograms
    std::map<std::string, ROOT::TThreadedObject<TH2D>> hsgas, hstwo, hszero;
    // Histogram models
    auto hGasSil{new TH2D{"hGasSil", ";E_{Sil} [MeV];#Delta E_{gas} [arb. units]", 450, 0, 70, 600, 0, 3000}};
    auto hTwoSils{new TH2D{"hTwoSils", ";#DeltaE_{0} [MeV];#DeltaE_{1} [MeV]", 500, 0, 80, 400, 0, 30}};
    for (const auto &layer : {"f0", "l0", "r0"})
    {
        hsgas.emplace(layer, *hGasSil);
        hsgas[layer]->SetTitle(TString::Format("%s", layer));
    }
    hstwo.emplace("f0-f1", *hTwoSils);
    hstwo["f0-f1"]->SetTitle("f0-f1");
    for (int s = 0; s < 4; s++)
    {
        hszero.emplace(std::to_string(s), *hTwoSils);
        hszero[std::to_string(s)]->SetTitle(TString::Format("f2_%d vs f3;E_{f3} [MeV];#DeltaE_{f2} [MeV]", s));
    }

    // Fill them
    df.Foreach(
        [&](ActRoot::MergerData &m)
        {
            // Light
            if (lambdaOne(m)) // Gas-E0 PID
            {
                auto layer{m.fLight.GetLayer(0)};
                if (hsgas.count(layer))
                    hsgas[layer]->Fill(m.fLight.fEs.front(), m.fLight.fQave);
            }
            else if (lambdaTwo(m)) // E0-E1 PID
            {
                hstwo["f0-f1"]->Fill(m.fLight.fEs[0], m.fLight.fEs[1]);
            }
            // Heavy
            if (m.fHeavy.GetNLayers() == 2)
            {
                auto n{std::to_string(m.fHeavy.fNs[0])};
                if (hszero.count(n))
                    hszero[n]->Fill(m.fHeavy.fEs[1], m.fHeavy.fEs[0]);
            }
        },
        {"MergerData"});

    // If cuts are present, apply them
    ActRoot::CutsManager<std::string> cuts;
    // Gas PID
    cuts.ReadCut("l0", TString::Format("./Cuts/pid_%s_l0.root", light.c_str()).Data());
    cuts.ReadCut("r0", TString::Format("./Cuts/pid_%s_r0.root", light.c_str()).Data());
    // cuts.ReadCut("f0", TString::Format("./Cuts/pid_%s_f0.root", light.c_str()).Data());
    // Two sils PID
    // cuts.ReadCut("f0-f1", TString::Format("./Cuts/pid_%s_f0_f1.root", light.c_str()).Data());
    if (cuts.GetListOfKeys().size())
    {
        // Apply PID and save in file
        auto gated{df.Filter(
            [&](ActRoot::MergerData &m)
            {
                // One silicon
                if (lambdaOne(m))
                {
                    auto layer{m.fLight.GetLayer(0)};
                    if (cuts.GetCut(layer))
                        return cuts.IsInside(layer, m.fLight.fEs[0], m.fLight.fQave);
                    else
                        return false;
                }
                else if (cuts.GetCut("f0-f1") && lambdaTwo(m))
                    return cuts.IsInside("f0-f1", m.fLight.fEs[0], m.fLight.fEs[1]);
                else
                    return false;
            },
            {"MergerData"})};
        auto name{TString::Format("./Outputs/tree_pid_%s_%s.root", target.c_str(), light.c_str())};
        std::cout << "Saving PID_Tree in file : " << name << '\n';
        gated.Snapshot("PID_Tree", name.Data());
    }

    // Draw
    auto *c0{new TCanvas{"c0", "PID canvas"}};
    c0->DivideSquare(6);
    int p{1};
    c0->cd(1);
    for (auto &[layer, h] : hsgas)
    {
        c0->cd(p);
        h.Merge()->DrawClone("colz");
        cuts.DrawCut(layer);
        p++;
    }
    for (auto &[layer, h] : hstwo)
    {
        c0->cd(p);
        h.Merge()->DrawClone("colz");
        p++;
    }

    auto *c1{new TCanvas{"c1", "PID canvas for 0 deg"}};
    c1->DivideSquare(hszero.size());
    p = 1;
    for (auto &[s, h] : hszero)
    {
        c1->cd(p);
        h.Merge()->DrawClone("colz");
        p++;
    }
}
