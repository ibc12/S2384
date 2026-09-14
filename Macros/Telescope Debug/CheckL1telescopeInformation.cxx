#include "ActCluster.h"
#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActLine.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActSRIM.h"
#include "ActSilData.h"
#include "ActSilSpecs.h"
#include "ActTPCData.h"
#include "ActTPCParameters.h"
#include "ActVoxel.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"
#include <random>

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMath.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TSpline.h>
#include <TStyle.h>

#include <Math/Point3D.h>
#include <Math/Vector3D.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

#include "../../PrettyStyle.C"


void CheckL1telescopeInformation()
{
    ROOT::EnableImplicitMT();

    ROOT::RDataFrame df_total {"PreProcessed_Tree", "../../PostAnalysis/Outputs/tree_preprocess_11Li.root"};
    ROOT::RDataFrame df {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_F_11Li_d_d_filtered.root"};
    auto def_L1 {df.Filter([](ActRoot::ModularData& m) { return m.Get("GATCONF") == 8; }, {"ModularData"})}; // only L1
    auto def_silicon =
        df.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == true; }, {"MergerData"}); // only silicon

    auto def_L1_total {df_total.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; },
                                       {"MergerData"})}; // only L1

    std::shared_ptr<ActPhysics::SilSpecs> sils = std::make_shared<ActPhysics::SilSpecs>();
    sils->ReadFile("../../configs/silspecs.conf");

    auto def_telescope = def_L1.Filter(
        [&sils](ActRoot::MergerData& m, ActRoot::SilData& s)
        {
            // Find if the map haskeys f2 and f3
            s.ApplyFinerThresholds(sils);

            auto it_f2 = s.fSiE.find("f2");
            auto it_f3 = s.fSiE.find("f3");

            bool has_f2 = it_f2 != s.fSiE.end() && !it_f2->second.empty();
            bool has_f3 = it_f3 != s.fSiE.end() && !it_f3->second.empty();

            return has_f2 && has_f3;
            // auto layers = m.fSilLayers;
            // for(auto layer : layers)
            // {
            //     std::cout << "Layer: " << layer << std::endl;
            // }
            // if(layers.size() != 2)
            //     return false;
            // for(auto layer : layers)
            // {
            //     if(layer != "f2" && layer != "f3")
            //     {
            //         return false;
            //     }
            // }
            // return layers.size() > 0;
        },
        {"MergerData", "SilData"}); // telescope for L1

    auto def_telescope_silicon_silData = def_silicon.Filter(
        [&sils](ActRoot::MergerData& m, ActRoot::SilData& s)
        {
            // Find if the map haskeys f2 and f3
            s.ApplyFinerThresholds(sils);

            auto it_f2 = s.fSiE.find("f2");
            auto it_f3 = s.fSiE.find("f3");

            bool has_f2 = it_f2 != s.fSiE.end() && !it_f2->second.empty();
            bool has_f3 = it_f3 != s.fSiE.end() && !it_f3->second.empty();

            return has_f2 && has_f3;
            // auto layers = m.fSilLayers;
            // if(layers.size() != 2)
            //     return false;
            // for(auto layer : layers)
            // {
            //     if(layer != "f2" && layer != "f3")
            //     {
            //         return false;
            //     }
            // }
            // return true;
        },
        {"MergerData", "SilData"}); // telescope for L1 with silicon

    auto def_telescope_silicon_mergerData = def_silicon.Filter(
        [](ActRoot::MergerData& m, ActRoot::SilData& s)
        {
            auto layers = m.fHeavy.fLayers;
            if(layers.size() != 2)
                return false;
            bool has_f2 = false;
            bool has_f3 = false;
            for(auto layer : layers)
            {
                if(layer == "f2")
                {
                    has_f2 = true;
                }
                if(layer == "f3")
                {
                    has_f3 = true;
                }
            }
            return has_f2 && has_f3;
        },
        {"MergerData", "SilData"}); // telescope for L1 with silicon

    auto def_telescope_silicon_SilData_noMergerData = def_telescope_silicon_silData.Filter(
        [](ActRoot::MergerData& m, ActRoot::SilData& s)
        {
            auto layers = m.fHeavy.fLayers;
            // if(layers.size() != 2)
            //     return false;
            // std::cout << "------------------------------" << std::endl;
            // std::cout << "Layers size: " << layers.size() << std::endl;
            bool has_f2 = false;
            bool has_f3 = false;
            for(auto layer : layers)
            {
                // std::cout << "Layer: " << layer << std::endl;
                if(layer == "f2")
                {
                    has_f2 = true;
                }
                if(layer == "f3")
                {
                    has_f3 = true;
                }
            }
            if(!(has_f2 && has_f3))
            {
                std::cout << "-------------------------------" << std::endl;
                auto E_f2 = s.fSiE["f2"].front();
                auto E_f3 = s.fSiE["f3"].front();
                std::cout << "Sildata E for f2: " << E_f2 << std::endl;
                std::cout << "Sildata E for f3: " << E_f3 << std::endl;
                std::cout << "-------------------------------" << std::endl;
            }
            // std::cout << "------------------------------" << std::endl;
            return !(has_f2 && has_f3);
        },
        {"MergerData", "SilData"}); // telescope for L1 with silicon

    auto def_telescope_total = def_L1_total.Filter(
        [&sils](ActRoot::MergerData& m, ActRoot::SilData& s)
        {
            // Find if the map haskeys f2 and f3
            s.ApplyFinerThresholds(sils);

            auto it_f2 = s.fSiE.find("f2");
            auto it_f3 = s.fSiE.find("f3");

            bool has_f2 = it_f2 != s.fSiE.end() && !it_f2->second.empty();
            bool has_f3 = it_f3 != s.fSiE.end() && !it_f3->second.empty();

            return has_f2 && has_f3;
            // auto layers = m.fSilLayers;
            // if(layers.size() != 2)
            //     return false;
            // for(auto layer : layers)
            // {
            //     if(layer != "f2" && layer != "f3")
            //     {
            //         return false;
            //     }
            // }
            // return true;
        },
        {"MergerData", "SilData"}); // telescope for L1 total

    // Save events that go through SilData but not MergerData
    // std::ofstream out("./Outputs/telescope_L1.dat");
    // def_telescope.Foreach([&](ActRoot::MergerData& m) { m.Stream(out); }, {"MergerData"});
    // out.close();

    std::cout << "Total events L1 deuterium: " << def_L1.Count().GetValue() << std::endl;
    // std::cout << "Total events L1 total: " << def_L1_total.Count().GetValue() << std::endl;
    std::cout << "Total events L1 deuterium with telescope: " << def_telescope.Count().GetValue() << std::endl;
    // std::cout << "Total events L1 total with telescope: " << def_telescope_total.Count().GetValue() << std::endl;
    std::cout << "Total events deuterium with silicon: " << def_silicon.Count().GetValue() << std::endl;
    std::cout << "Total events deuterium with telescope and silicon (SilData): "
              << def_telescope_silicon_silData.Count().GetValue() << std::endl;
    std::cout << "Total events deuterium with telescope and silicon (MergerData): "
              << def_telescope_silicon_mergerData.Count().GetValue() << std::endl;
    std::cout << "Total events deuterium with telescope and silicon (SilData no MergerData): "
              << def_telescope_silicon_SilData_noMergerData.Count().GetValue() << std::endl;

    // Plot the heavy PID for the telescope L1 events
    // First define f2Energy and f3Energy in the df
    auto df_telescope_with_energy = def_telescope
                                        .Define("f2Energy",
                                                [](ActRoot::SilData& s)
                                                {
                                                    auto it_f2 = s.fSiE.find("f2");
                                                    if(it_f2 != s.fSiE.end() && !it_f2->second.empty())
                                                        return it_f2->second.front();
                                                    else
                                                        return -1.0f;
                                                },
                                                {"SilData"})
                                        .Define("f3Energy",
                                                [](ActRoot::SilData& s)
                                                {
                                                    auto it_f3 = s.fSiE.find("f3");
                                                    if(it_f3 != s.fSiE.end() && !it_f3->second.empty())
                                                        return it_f3->second.front();
                                                    else
                                                        return -1.0f;
                                                },
                                                {"SilData"});

    // Now plot f2 in y axis and f3 in x axis
    auto h_f2_f3 = df_telescope_with_energy.Histo2D(
        {"h_f2_f3", "f2 vs f3 energy; f3 energy [MeV]; f2 energy [MeV]", 300, 0, 80, 100, 0, 20}, "f3Energy",
        "f2Energy");
    ActRoot::CutsManager<std::string> cuts;
    cuts.ReadCut("f2_f3", "../../PostAnalysis/Cuts/pid_11Li_f2_0_11Li.root");
    auto c_f2_f3 = new TCanvas("c_f2_f3", "f2 vs f3 energy", 800, 600);
    h_f2_f3->DrawClone("colz");
    cuts.DrawCut("f2_f3");

    // Save events in the inverse banana of the heavy PID
    // std::ofstream out("./Outputs/telescope_L1_PID_heavy_inverse_banana.dat");
    // df_telescope_with_energy.Foreach(
    //     [&](ActRoot::MergerData& m, float f3Energy, float f2Energy)
    //     {
    //         if(f3Energy < 30 && f3Energy > 0 && f2Energy > 0 && f2Energy < 6)
    //             m.Stream(out);
    //     },
    //     {"MergerData", "f3Energy", "f2Energy"});
    // out.close();
    //
    // std::ofstream out1("./Outputs/telescope_L1_PID_heavy_beam.dat");
    // df_telescope_with_energy.Foreach(
    //     [&](ActRoot::MergerData& m, float f3Energy, float f2Energy)
    //     {
    //         if(f3Energy < 72 && f3Energy > 68 && f2Energy > 6 && f2Energy < 7.5)
    //             m.Stream(out1);
    //     },
    //     {"MergerData", "f3Energy", "f2Energy"});
    // out1.close();
}