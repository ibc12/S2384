
#include "ActModularData.h"
#include "ActSilData.h"

#include "RtypesCore.h"
#include <ROOT/RDataFrame.hxx>

#include "TCanvas.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"

#include <algorithm>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <vector>

void checkTimeStampSIlicons()
{
    // Get data from file in TTree
    auto fileSilicon {new TFile("../RootFiles/Data/Data_Run_0004.root")};
    auto treeSilicon {fileSilicon->Get<TTree>("VXITree")};

    ROOT::RDataFrame df(*treeSilicon);


    auto def = df.Define("Energy",
                         [&](ActRoot::SilData& d)
                         {
                             std::string layer {"l0"};
                             if(d.fSiN.count(layer))
                             {
                                 auto& ns {d.fSiN[layer]};
                                 auto it {std::find(ns.begin(), ns.end(), 0)};
                                 if(it != ns.end())
                                 {
                                     auto index = std::distance(ns.begin(), it);
                                     return d.fSiE[layer][index];
                                 }
                             }

                             return -1.f;
                         },
                         {"SilData"})
                   .Define("time",
                           [&](ActRoot::ModularData& md)
                           {
                               auto timeh_up {md.Get("CTR_TIMEH_UP")};
                               auto timeh {md.Get("CTR_TIMEH")};
                               auto timeml_up {md.Get("CTR_TIMEML_UP")};
                               auto timeml {md.Get("CTR_TIMEML")};

                               return timeml + timeml_up * 2e16 + timeh * 2e32 + timeh_up * 2e48;
                           },
                           {"ModularData"})
                   .Alias("timeCorrected", "time"); // Workaround until we dont have correct TSs

    // auto minTime {*def.Min("time")};
    // def = def.Define("timeCorrected",
    //                  [&](double& time)
    //                  {
    //                      // Correct the time
    //                      return (time - minTime);
    //                  },
    //                  {"time"});

    // Create the th2 model
    ROOT::RDF::TH2DModel mTSE {"hTSE", "TS vs E;TS [10 #mus]; E_{Sil} [channel]", 600, 0, 600e7, 1000, 0, 5000};
    ROOT::RDF::TH2DModel mEntryE {"hEntryE", "Entry vs E;Entry; E_{Sil} [channel]", 600, 0, 2e6, 16384, 0, 16384};
    auto hTSE {def.Histo2D(mTSE, "timeCorrected", "Energy")};
    auto hEntryE {def.Histo2D(mEntryE, "rdfentry_", "Energy")};
    auto hTS {def.Histo1D("timeCorrected")};
    auto hml {def.Define("ml", "fLeaves[\"CTR_TIMEML\"]").Histo1D("ml")};

    // Book histograms per layer and silicons
    std::map<std::string, std::vector<std::shared_ptr<TH2D>>> hs;
    int nsil {12};
    for(const auto& layer : {"l0", "r0", "f0"})
    {
        hs[layer] = {};
        for(int s = 0; s < nsil; s++)
        {
            auto h {mEntryE.GetHistogram()};
            h->SetDirectory(nullptr);
            auto aux {TString(layer)};
            aux.ToUpper();
            auto name {aux + "_" + TString::Format("%d", s)};
            h->SetNameTitle(name, name);
            hs[layer].push_back(h);
        }
    }
    // Fill!
    def.Filter("rdfentry_ < 230000")
        .Foreach(
            [&](ULong64_t entry, ActRoot::SilData& d)
            {
                for(const auto& layer : d.GetLayers())
                {
                    if(!hs.count(layer))
                        continue;
                    auto mult {d.GetMult(layer)};
                    for(int i = 0; i < mult; i++)
                    {
                        auto n {d.fSiN[layer][i]};
                        auto e {d.fSiE[layer][i]};
                        // std::cout << layer << " n : " << n << " e : " << e << '\n';
                        hs[layer][n]->Fill(entry, e);
                    }
                }
            },
            {"rdfentry_", "SilData"});

    // Write them
    auto fout {std::make_shared<TFile>("./Outputs/SiWall_gated_entry_run4.root", "recreate")};
    // F0
    fout->mkdir("Raw/F0");
    fout->cd("Raw/F0");
    for(const auto& h : hs["f0"])
    {
        auto* py {h->ProjectionY(h->GetName())};
        py->Write();
        delete py;
    }
    // L0
    fout->cd();
    fout->mkdir("Raw/L0");
    fout->cd("Raw/L0");
    for(const auto& h : hs["l0"])
    {
        auto* py {h->ProjectionY(h->GetName())};
        py->Write();
        delete py;
    }
    // R0
    fout->cd();
    fout->mkdir("Raw/R0");
    fout->cd("Raw/R0");
    for(const auto& h : hs["r0"])
    {
        auto* py {h->ProjectionY(h->GetName())};
        py->Write();
        delete py;
    }
    fout->Close();

    auto c0 {new TCanvas("c0", "Silicon Energy vs time")};
    c0->DivideSquare(2);
    c0->cd(1);
    hTSE->DrawClone("colz");
    c0->cd(2);
    hEntryE->DrawClone("colz");

    // auto cTime {new TCanvas("cTime", "Silicon Time")};
    // cTime->DivideSquare(4);
    // cTime->cd(1);
    // hTS->DrawClone("colz");
    // cTime->cd(2);
    // hml->DrawClone();
    for(const auto& [layer, vh] : hs)
    {
        auto c {new TCanvas {TString::Format("c%s", layer.c_str()), TString::Format("Canvas for %s", layer.c_str())}};
        c->DivideSquare(vh.size());
        for(int i = 0; i < vh.size(); i++)
        {
            c->cd(i + 1);
            vh[i]->DrawClone("colz");
        }
    }
}
