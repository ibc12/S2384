#include "TCanvas.h"
#include "TColor.h"

#include "ActMergerData.h"
#include "ActCutsManager.h"

#include "ROOT/RDataFrame.hxx"

void GetHeavyGates()
{
    ROOT::DisableImplicitMT();
    ROOT::RDataFrame df {"Final_Tree", "../PostAnalysis/Outputs/tree_ex_11Li_d_d_filtered.root"};


    std::string beam {"11Li"};
    std::vector<std::string> listOfCuts {"9Li", "11Li"};

    // 1st condition -> Impose light hits the silicon (otherwise L1 trigger doesnt have Heavy hit either)
    // 2nd condition -> Ensure f2 AND f3
    auto gated {df.Filter([](ActRoot::MergerData& mer)
                          { return mer.fLight.IsFilled() && mer.fHeavy.GetNLayers() == 2; }, {"MergerData"})};


    // Read cuts files for heavy particle, there is one per pad of f2 telescope (divide in 4 squares)
    ActRoot::CutsManager<std::string> cuts;
    // Read cuts on heavy particle
    for(const auto& recoil : listOfCuts)
    {
        for(int s {0}; s < 4; s++)
        {
            auto cutfile {TString::Format("../PostAnalysis/Cuts/pid_%s_f2_%d_11Li.root", recoil.c_str(), s)};
            cuts.ReadCut(TString::Format("%s_%d", recoil.c_str(), s).Data(), cutfile.Data());
        }
    }

    // Filter de original dataframe with each cut checking whether the particle is inside it

    auto df9Li = gated.Filter(
        [&cuts](ActRoot::MergerData& mer)
        {
            auto n {std::to_string(mer.fHeavy.fNs[0])};
            auto key {TString::Format("9Li_%s", n.c_str())};
            return cuts.IsInside(key.Data(), mer.fHeavy.fEs[1], mer.fHeavy.fEs[0]);
        },
        {"MergerData"});

    auto df11Li = gated.Filter(
        [&cuts](ActRoot::MergerData& mer)
        {
            auto n {std::to_string(mer.fHeavy.fNs[0])};
            auto key {TString::Format("11Li_%s", n.c_str())};
            return cuts.IsInside(key.Data(), mer.fHeavy.fEs[1], mer.fHeavy.fEs[0]);
        },
        {"MergerData"});

    // Save the filtered dataframes as you want. Here I just plot it to  show that it works
    auto c {new TCanvas("cHeavyGates", "cHeavyGates", 800, 600)};
    c->Divide(3,1);
    c->cd(1);
    auto hExAll = gated.Histo1D({"hExAfter", "Excitation Energy after filtering;E_{x} (MeV);Counts", 100, -5, 10}, "Ex");
    hExAll->DrawClone();
    c->cd(2);
    auto hEx9Li = df9Li.Histo1D({"hEx9Li", "Excitation Energy for 9Li;E_{x} (MeV);Counts", 100, -5, 10}, "Ex");
    hEx9Li->DrawClone();
    c->cd(3);
    auto hEx11Li = df11Li.Histo1D({"hEx11Li", "Excitation Energy for 11Li;E_{x} (MeV);Counts", 100, -5, 10}, "Ex");
    hEx11Li->DrawClone();

        
    // Savein .root only the branch Ex and ThetaCM for each filtered dataframe
    df9Li.Snapshot("HeavyGated9Li", "./Outputs/tree_heavy_9Li_11Li_d_filtered.root", {"Ex", "ThetaCM"});
    df11Li.Snapshot("HeavyGated11Li", "./Outputs/tree_heavy_11Li_11Li_d_filtered.root", {"Ex", "ThetaCM"});


    
}