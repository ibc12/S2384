#include "TCanvas.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TGraph2DAsymmErrors.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TTree.h"

#include <iostream>
#include <string>
#include <vector>

void CompareEfficienciesDisplacementBEamTelescope()
{
    std::vector<std::string> files = {
        "../../Simulation/Outputs/11Li/test_telescope_beam_displacement/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_1.0.root",
        "../../Simulation/Outputs/11Li/test_telescope_beam_displacement/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_3.0.root",
        "../../Simulation/Outputs/11Li/test_telescope_beam_displacement/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_7.0.root",
        "../../Simulation/Outputs/11Li/test_telescope_beam_displacement/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_15.0.root",
        "../../Simulation/Outputs/11Li/test_telescope_beam_displacement/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_20.0.root",
        "../../Simulation/Outputs/11Li/test_telescope_beam_displacement/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_25.0.root",
        "../../Simulation/Outputs/11Li/test_telescope_beam_displacement/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_30.0.root"
    };

    std::vector<TEfficiency*> efficiencies;
    std::vector<TFile*> openedFiles;

    for(const auto& file : files)
    {
        TFile* f = TFile::Open(file.c_str(), "READ");
        if(!f || f->IsZombie())
        {
            std::cout << "Cannot open " << file << std::endl;
            continue;
        }

        auto* eff = dynamic_cast<TEfficiency*>(f->Get("effCMgateHeavy"));
        if(!eff)
        {
            std::cout << "Cannot find effCM in " << file << std::endl;
            continue;
        }

        efficiencies.push_back(eff);
        openedFiles.push_back(f); // mantener el fichero abierto
    }

    TCanvas* c = new TCanvas("c", "Efficiencies", 1000, 700);

    std::vector<int> colors = {kBlack,      kRed,      kBlue,   kGreen + 2, kMagenta,
                               kOrange + 1, kCyan + 2, kViolet, kPink + 7,  kGray + 2};

    TLegend* legend = new TLegend(0.60, 0.15, 0.88, 0.45);
    legend->SetBorderSize(0);

    for(size_t i = 0; i < efficiencies.size(); ++i)
    {
        efficiencies[i]->SetLineColor(colors[i % colors.size()]);
        efficiencies[i]->SetMarkerColor(colors[i % colors.size()]);
        efficiencies[i]->SetMarkerStyle(20 + i);

        // Get the labels by the last word between .root and the last _
        TString label;
        TString fileName = files[i];
        Ssiz_t posRoot = fileName.Last('.');
        Ssiz_t posUnderscore = fileName.Last('_');
        if(posUnderscore != kNPOS && posRoot != kNPOS)
        {
            label = fileName(posUnderscore + 1, posRoot - posUnderscore - 1);
        }

        if(i == 0)
        {
            efficiencies[i]->Draw("AP");

            gPad->Update();

            auto* h = efficiencies[i]->GetPaintedHistogram();
            if(h)
            {
                h->GetXaxis()->SetRangeUser(0, 30);
                h->GetYaxis()->SetRangeUser(0, 1.05);
                h->GetXaxis()->SetTitle("#theta_{CM} [deg]");
                h->GetYaxis()->SetTitle("Efficiency");
            }
        }
        else
        {
            efficiencies[i]->Draw("P SAME");
        }

        legend->AddEntry(efficiencies[i], label, "lp");
    }

    legend->Draw();

    c->Modified();
    c->Update();
}