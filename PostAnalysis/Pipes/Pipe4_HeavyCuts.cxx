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

#include "../../PrettyStyle.C"
#include "../HistConfig.h"

void Pipe4_HeavyCuts(const std::string& beam, const std::string& target, const std::string& light)
{
    PrettyStyle(false);
    // gStyle->SetPalette(kViridis);

    auto infile {
        TString::Format("./Outputs/tree_ex_%s_%s_%s_filtered.root", beam.c_str(), target.c_str(), light.c_str())};
    // ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"Final_Tree", infile.Data()};

    std::vector<std::string> listOfCuts {"7Li"};

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
    for(int i = 1; i < hsEx.size(); ++i) // empezamos en 1 para saltar "All"
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

    // Plots Ex tot vs Ex gated for each cut

    std::vector<int> cutColors = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 1};

    for(size_t i = 1; i < hsEx.size(); ++i) // empezamos en 1 (saltamos "All")
    {
        const auto& cutName = labels[i];

        auto* c =
            new TCanvas(TString::Format("cTot_%s", cutName.c_str()), TString::Format("Total vs %s", cutName.c_str()));
        c->cd();

        // Histograma total
        auto* hTot = (TH1D*)hsEx[0].GetPtr()->Clone(TString::Format("hTot_%s", cutName.c_str()));
        hTot->SetLineColor(kBlack);
        hTot->SetLineWidth(2);

        // Histograma del corte
        auto* hCut = (TH1D*)hsEx[i].GetPtr()->Clone(TString::Format("h_%s", cutName.c_str()));
        hCut->SetLineColor(cutColors[(i - 1) % cutColors.size()]);
        hCut->SetLineWidth(2);

        // Dibujar
        hTot->Draw("hist");
        hCut->Draw("hist same");

        // Leyenda
        auto* leg = new TLegend(0.6, 0.7, 0.88, 0.88);
        leg->AddEntry(hTot, "Total", "l");
        leg->AddEntry(hCut, cutName.c_str(), "l");
        leg->Draw();
    }

    // Kinematic plot
    auto* cKin {new TCanvas {"c42", "Pipe4 Kinematics"}};
    cKin->cd();
    auto hKinClone = (TH2D*)hKin.GetPtr()->Clone("hKinClone");
    hKinClone->Draw("colz");
}

#endif
