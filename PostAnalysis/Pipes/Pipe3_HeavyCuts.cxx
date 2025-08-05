#ifndef PIPE3_HEAVYCUTS_H
#define PIPE3_HEAVYCUTS_H
#include "ActCutsManager.h"
#include "ActMergerData.h"

#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RResultPtr.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TColor.h"
#include "THStack.h"
#include "TString.h"
#include "TStyle.h"
#include "TVirtualPad.h"

#include <string>
#include <vector>

#include "../HistConfig.h"

void Pipe3_HeavyCuts(const std::string& beam, const std::string& target, const std::string& light)
{
    gStyle->SetPalette(kViridis);

    auto infile {TString::Format("./Outputs/tree_ex_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str())};
    // ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"Final_Tree", infile.Data()};

    // Apply cuts on
    // 1-> Impose light hits the silicon (otherwise L1 trigger doesnt have Heavy hit either)
    // 2-> Ensure f2 AND f3
    auto gated {df.Filter([](ActRoot::MergerData& mer)
                          { return mer.fLight.IsFilled() && mer.fHeavy.GetNLayers() == 2; }, {"MergerData"})};

    // Read cuts for heavy particle
    ActRoot::CutsManager<std::string> cuts;
    // Read cuts on heavy particle
    for(const auto& recoil : {"9Li", "11Li"})
    {
        for(int s {0}; s < 4; s++)
        {
            auto cutfile {TString::Format("./Cuts/pid_%s_f2_%d_%s.root", recoil, s, beam.c_str())};
            cuts.ReadCut(TString::Format("%s_%d", recoil, s).Data(), cutfile.Data());
        }
    }

    // Declare lambda
    // Apply gates
    std::vector<ROOT::RDF::RNode> nodes {gated};
    for(const auto& particle : {"9Li", "11Li"})
    {
        std::string selection {particle};
        nodes.push_back(gated.Filter(
            [selection, &cuts](ActRoot::MergerData& mer)
            {
                auto key {selection + "_" + std::to_string(mer.fHeavy.fNs[0])};
                if(cuts.IsInside(key, mer.fHeavy.fEs[1], mer.fHeavy.fEs[0]))
                    return true;
                else
                    return false;
            },
            {"MergerData"}));
    }
    // And book histograms
    std::vector<std::string> labels {"All", "^{9}Li", "^{11}Li"};
    std::vector<ROOT::RDF::RResultPtr<TH1D>> hsEx;
    int idx {};
    for(auto& node : nodes)
    {
        auto hEx {node.Histo1D(HistConfig::Ex, "Ex")};
        hEx->SetTitle(labels[idx].c_str());

        hsEx.push_back(hEx);
        idx++;
    }

    std::map<std::string, ROOT::TThreadedObject<TH2D>> hsgas, hstwo, hszero;
    auto hTwoSils {new TH2D {"hTwoSils", ";#DeltaE_{0} [MeV];#DeltaE_{1} [MeV]", 500, 0, 80, 400, 0, 30}};
    for(int s = 0; s < 4; s++)
    {
        hszero.emplace(std::to_string(s), *hTwoSils);
        hszero[std::to_string(s)]->SetTitle(TString::Format("f2_%d vs f3;E_{f3} [MeV];#DeltaE_{f2} [MeV]", s));
    }
    gated.Foreach(
        [&](ActRoot::MergerData& mer)
        {
            auto n {std::to_string(mer.fHeavy.fNs[0])};
            if(hszero.count(n))
                hszero[n]->Fill(mer.fHeavy.fEs[1], mer.fHeavy.fEs[0]);
        },
        {"MergerData"});


    // Create stack
    auto* stack {new THStack};
    stack->SetTitle(HistConfig::Ex.fTitle);
    for(auto& h : hsEx)
    {
        h->SetLineWidth(1);
        stack->Add((TH1D*)h->Clone());
    }

    // Draw
    auto* c1 {new TCanvas {"c1", "Pipe3 Heavy PID"}};
    c1->DivideSquare(hszero.size());
    int p {1};
    for(auto& [s, h] : hszero)
    {
        c1->cd(p);
        h.Merge()->DrawClone("colz");
        for(const auto& particle : {"11Li", "9Li"})
        {
            auto key {TString::Format("%s_%s", particle, s.c_str())};
            cuts.DrawCut(key.Data());
        }
        p++;
    }

    auto* c0 {new TCanvas {"c0", "Pipe3 Ex gated"}};
    // c0->DivideSquare(4);
    // c0->cd(1);
    stack->Draw("nostack plc");
    stack->GetXaxis()->SetRangeUser(-5, 12);

    gPad->BuildLegend();
}

#endif
