#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "THStack.h"
#include "TString.h"
#include "TStyle.h"
#include "TVirtualPad.h"

#include <memory>
#include <string>
#include <vector>
void compOldNew()
{
    gStyle->SetPalette(kViridis);

    std::vector<std::string> channels {"p", "d", "t"};
    std::vector<THStack*> hsSil, hsL1;

    for(const auto& channel : channels)
    {
        auto path {TString::Format("../PostAnalysis/Outputs/")};
        hsSil.push_back(new THStack);
        hsSil.back()->SetTitle(TString::Format("^{11}Li(d,%s) silicon;E_{x} [MeV];Counts / 150 keV", channel.c_str()));
        hsL1.push_back(new THStack);
        hsL1.back()->SetTitle(TString::Format("^{11}Li(d,%s) L1;E_{x} [MeV];Counts / 150 keV", channel.c_str()));
        int idx {};
        for(const auto& subdir : {"BeginExp", ""})
        {
            auto file {std::make_shared<TFile>(
                TString::Format("%s/%s/tree_ex_11Li_d_%s.root", path.Data(), subdir, channel.c_str()).Data())};
            auto hSil {file->Get<TH1D>("hExSil")};
            auto hL1 {file->Get<TH1D>("hExL1")};
            hSil->SetDirectory(nullptr);
            hL1->SetDirectory(nullptr);
            for(auto* h : {hSil, hL1})
            {
                h->SetTitle(idx == 0 ? "Before" : "After");
                h->Rebin(2);
                h->SetLineWidth(2);
            }

            // Push to stack
            hsSil.back()->Add(hSil);
            hsL1.back()->Add(hL1);
            idx++;
        }
    }
    // Sum
    for(auto* vec : {&hsSil, &hsL1})
    {
        for(auto* stack : *vec)
        {
            int idx {};
            TH1D* clone {};
            for(int i = 0; i < stack->GetNhists(); i++)
            {
                if(i == 0)
                    clone = (TH1D*)stack->GetHists()->At(i)->Clone();
                else
                    clone->Add((TH1D*)stack->GetHists()->At(i));
            }
            clone->SetTitle("Sum");
            stack->Add(clone);
        }
    }

    // Draw
    auto* c0 {new TCanvas {"c0", "Sil canvas"}};
    c0->DivideSquare(hsSil.size());
    for(int i = 0; i < hsSil.size(); i++)
    {
        c0->cd(i + 1);
        hsSil[i]->Draw("nostack plc");
        hsSil[i]->GetXaxis()->SetRangeUser(-2, 12);
        gPad->BuildLegend();
    }

    auto* c1 {new TCanvas {"c1", "L1 canvas"}};
    c1->DivideSquare(hsL1.size());
    for(int i = 0; i < hsL1.size(); i++)
    {
        c1->cd(i + 1);
        hsL1[i]->Draw("nostack plc");
        hsL1[i]->GetXaxis()->SetRangeUser(-2, 12);
        gPad->BuildLegend();
    }
}
