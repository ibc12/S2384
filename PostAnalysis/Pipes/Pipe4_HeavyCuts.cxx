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

#include "../../PrettyStyle.C"

void Pipe4_HeavyCuts(const std::string& beam, const std::string& target, const std::string& light)
{
    PrettyStyle(false);
    //gStyle->SetPalette(kViridis);

    auto infile {TString::Format("./Outputs/tree_ex_%s_%s_%s_filtered.root", beam.c_str(), target.c_str(), light.c_str())};
    // ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"Final_Tree", infile.Data()};

    std::vector<std::string> listOfCuts {"9Li", "11Li"};

    // Apply cuts on
    // 1-> Impose light hits the silicon (otherwise L1 trigger doesnt have Heavy hit either)
    // 2-> Ensure f2 AND f3
    auto gated {df.Filter([](ActRoot::MergerData& mer)
                          { return mer.fLight.IsFilled() && mer.fHeavy.GetNLayers() == 2; }, {"MergerData"})};

    // Read cuts for heavy particle
    ActRoot::CutsManager<std::string> cuts;
    // Read cuts on heavy particle
    for(const auto& recoil : listOfCuts)
    {
        for(int s {0}; s < 4; s++)
        {
            auto cutfile {TString::Format("./Cuts/pid_%s_f2_%d_%s.root", recoil.c_str(), s, beam.c_str())};
            cuts.ReadCut(TString::Format("%s_%d", recoil.c_str(), s).Data(), cutfile.Data());
        }
    }

    // Declare lambda
    // Apply gates
    std::vector<ROOT::RDF::RNode> nodes {gated};
    for(const auto& particle : listOfCuts)
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
    std::vector<std::string> labels {"All"};
    for(const auto& cut : listOfCuts)
        labels.push_back(cut);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> hsEx;
    int idx {};
    for(auto& node : nodes)
    {
        auto hEx {node.Histo1D(HistConfig::Ex, "Ex")};
        hEx->SetTitle(labels[idx].c_str());

        hsEx.push_back(hEx);
        idx++;
    }

    // Do a rebining for a factor 2, to see better the results
    for(auto& h : hsEx)
    {
        h->Rebin(2); // rebin factor 2

        // Actualizar título con el nuevo bin width
        double binWidth_keV = ((35. - (-10.)) / 600 * 2) * 1e3; // 2*original width en keV
        TString newTitle = TString::Format("Excitation energy;E_{x} [MeV];Counts / %.f keV", binWidth_keV);
        h->SetTitle(newTitle);
    }

    std::map<std::string, ROOT::TThreadedObject<TH2D>> hsgas, hstwo, hszero;
    auto hTwoSils_Pipe3 {new TH2D {"hTwoSils3", ";#DeltaE_{0} [MeV];#DeltaE_{1} [MeV]", 500, 0, 80, 400, 0, 30}};
    for(int s = 0; s < 4; s++)
    {
        hszero.emplace(std::to_string(s), *hTwoSils_Pipe3);
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
    double binWidth_keV = ((35. - (-10.)) / 600 * 2) * 1e3; // 2*original width en keV
    TString newTitle = TString::Format("Excitation energy;E_{x} [MeV];Counts / %.f keV", binWidth_keV);
    stack->SetTitle(newTitle);
    for(auto& h : hsEx)
    {
        h->SetLineWidth(1);
        stack->Add((TH1D*)h->Clone());
    }

    // Get kinematic plot on gated in heavy detection
    auto hKin {gated.Histo2D(HistConfig::KinEl, "fThetaLight", "EVertex")};

    // Draw
    auto* c1 {new TCanvas {"c31", "Pipe3 Heavy PID"}};
    c1->DivideSquare(hszero.size());
    int p {1};
    for(auto& [s, h] : hszero)
    {
        c1->cd(p);
        h.Merge()->DrawClone("colz");
        for(const auto& particle : listOfCuts)
        {
            auto key {TString::Format("%s_%s", particle.c_str(), s.c_str())};
            cuts.DrawCut(key.Data());
        }
        p++;
    }

    // Canvas for each Ex spectrum
    for(int i = 1; i < hsEx.size(); ++i)  // empezamos en 1 para saltar "All"
    {
        auto* c = new TCanvas(TString::Format("c3_%s", labels[i].c_str()),
                            TString::Format("Excitation energy %s", labels[i].c_str()));
        hsEx[i]->DrawClone("hist");
    }

    auto* c0 {new TCanvas {"c40", "Pipe4 Ex gated"}};
    // c0->DivideSquare(4);
    // c0->cd(1);
    stack->Draw("nostack plc");
    stack->GetXaxis()->SetRangeUser(-5, 12);

    gPad->BuildLegend();
    // ================
    // Canvas extra 1: Solo espectro total
    // ================
    auto* cTot = new TCanvas("cTot", "Espectro total");
    cTot->cd();

    auto* hTot = (TH1D*)hsEx[0].GetPtr()->Clone("hTot");
    hTot->SetLineWidth(2);
    hTot->SetLineColor(kBlack);
    hTot->Draw("hist");

    // ================
    // Canvas extra 2: Total + 9Li
    // ================
    auto* cTot_9Li = new TCanvas("cTot_9Li", "Total vs 9Li");
    cTot_9Li->cd();

    auto* hTot2 = (TH1D*)hsEx[0].GetPtr()->Clone("hTot2");
    hTot2->SetLineColor(kBlack);
    hTot2->SetLineWidth(2);

    auto* h9 = (TH1D*)hsEx[1].GetPtr()->Clone("h9_clone");
    h9->SetLineColor(kRed+1);
    h9->SetLineWidth(2);

    hTot2->Draw("hist");
    h9->Draw("hist same");

    auto* legend9 = new TLegend(0.6,0.7,0.88,0.88);
    legend9->AddEntry(hTot2, "Total", "l");
    legend9->AddEntry(h9, "9Li", "l");
    legend9->Draw();

    // ================
    // Canvas extra 3: Total + 11Li
    // ================
    auto* cTot_11Li = new TCanvas("cTot_11Li", "Total vs 11Li");
    cTot_11Li->cd();

    auto* hTot3 = (TH1D*)hsEx[0].GetPtr()->Clone("hTot3");
    hTot3->SetLineColor(kBlack);
    hTot3->SetLineWidth(2);

    auto* h11 = (TH1D*)hsEx[2].GetPtr()->Clone("h11_clone");
    h11->SetLineColor(kBlue+1);
    h11->SetLineWidth(2);

    hTot3->Draw("hist");
    h11->Draw("hist same");

    auto* legend11 = new TLegend(0.6,0.7,0.88,0.88);
    legend11->AddEntry(hTot3, "Total", "l");
    legend11->AddEntry(h11, "11Li", "l");
    legend11->Draw();

    // Kinematic plot
    auto* cKin {new TCanvas {"c42", "Pipe4 Kinematics"}};
    cKin->cd();
    auto hKinClone = (TH2D*)hKin.GetPtr()->Clone("hKinClone");
    hKinClone->Draw("colz");

}

#endif
