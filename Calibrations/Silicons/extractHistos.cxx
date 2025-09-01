#include "TChain.h"
#include "ROOT/RDataFrame.hxx"
#include "ActSilData.h"
#include "TCanvas.h"

void extractHistos(){

    auto* chain = new TChain("VXITree");
    chain->Add("../../RootFiles/Data/Data_Run_0130.root");
    ROOT::RDataFrame df {*chain};

    std::map<std::string, std::vector<TH1D*>> hs;
    int nsils {};
    for(const auto& layer : {"r0", "l0", "f0", "f2"})
    {
        if(layer == "f2")
            nsils = 4; // f2 has only 4 pads
        else
            nsils = 12;
        for(int i =0; i < nsils; i++){
            auto* h = new TH1D(TString::Format("h%s%d", layer, i), TString::Format("%s_%d;Channel;Counts", layer, i), 16384, 0, 16384);
            auto aux {TString(layer)};
            aux.ToUpper();
            auto name {aux+"_"+TString::Format("%d",i)};
            h->SetNameTitle(name,name);
            hs[layer].push_back(h);
        }
    }

    df.Foreach([&](ActRoot::SilData& sil){ // the [&] is to inherit the elemets from the previous for loop
        for(const auto& [layer, ns] : sil.fSiN){
            if(!hs.count(layer))
                continue;//account for f1 that doean't exist
            for(int i = 0; i < ns.size(); i++){
                auto n = ns[i];
                auto e = sil.fSiE[layer][i];
                hs[layer][i]->Fill(e);
            }
        }
    }, {"SilData"});

    //Draw one layer per canvas
    for(auto& [layer, vec] : hs){
        auto* c = new TCanvas(TString::Format("c%s", layer.c_str()), layer.c_str());
        c->DivideSquare(vec.size());
        for(int i =0;i <vec.size(); i++){
            c->cd(i + 1);
            vec[i]->Draw();
        }
    }

    //Write histograms to file
    auto fout {std::make_shared<TFile>("./Inputs/Si_calib_histos_run0130.root", "recreate")};
    // F0
    fout->mkdir("Raw/F0");
    fout->cd("Raw/F0");
    for(const auto& h : hs["f0"]){
        h->Write();
    }
    // L0
    fout->mkdir("Raw/L0");
    fout->cd("Raw/L0");
    for(const auto& h : hs["l0"]){
        h->Write();
    }
    // R0
    fout->mkdir("Raw/R0");
    fout->cd("Raw/R0");
    for(const auto& h : hs["r0"]){
        h->Write();
    }
    // F2
    fout->mkdir("Raw/F2");
    fout->cd("Raw/F2");
    for(const auto& h : hs["f2"]){
        h->Write();
    }

    fout->Close();
}