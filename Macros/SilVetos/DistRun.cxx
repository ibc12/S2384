#include "ActDataManager.h"
#include "ActLine.h"
#include "ActMergerData.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TFile.h"
#include "TH1.h"
#include "TString.h"

#include "Math/Point3Dfwd.h"

#include <memory>
#include <tuple>
#include <utility>
#include <vector>

// 1st macro to execute

void DistRun()
{
    ROOT::EnableImplicitMT();

    ActRoot::DataManager data {"../../configs/data.conf", ActRoot::ModeType::EMerge};
    // data.SetRuns(19, 66); // L0 trigger window changed for
    data.SetRuns(68, 122); // L0 trigger window changed for
    auto chain {data.GetJoinedData()};

    ROOT::RDataFrame df {*chain};

    std::string layer {"f0"};

    // Filter side events
    auto gated {df.Filter(
        [&](ActRoot::MergerData& merger)
        {
            if(merger.fSilLayers.size() == 1)
                if(merger.fSilLayers.front() == layer)
                    return true;
            return false;
        },
        {"MergerData"})};

    //  Define distances in mm
    double base {0.};
    std::vector<double> dists;
    // l0 from 310 to 320
    // f0 from 315 to 330
    for(double d = 305; d < 315; d += 1)
        dists.push_back(base + d);

    int xbins {200};
    int zbins {};
    if(layer == "f0")
        zbins = 250;
    else
        zbins = 200;
    std::pair<double, double> xlims {};
    std::pair<double, double> zlims {};
    if(layer == "l0")
    {
        xlims = {-20, 300};
        zlims = {130, 450};
    }
    else if(layer == "f0")
    {
        xlims = {-40, 290};
        zlims = {80, 500};
    }
    else if(layer == "r0")
    {
        xlims = {-20, 300};
        zlims = {150, 400};
    }


    // Save
    // TString outpath {TString::Format("./Outputs/Dists/histos_%s_preL0change.root", layer.c_str())};
    TString outpath {TString::Format("./Outputs/Dists/histos_%s.root", layer.c_str())};
    auto f {std::make_unique<TFile>(outpath, "recreate")};
    f->WriteObject(&dists, "dists");
    for(const auto& dist : dists)
    {
        std::cout << "Distance : " << dist << '\n';
        // Redefine
        auto node {gated.Define("NewSP",
                                [&, dist](ActRoot::MergerData& d)
                                {
                                    // Reconstruct line from BP and SP
                                    auto p {d.fBP};
                                    auto dir {(d.fSP - d.fBP)};
                                    ActRoot::Line line {p, dir, 0};
                                    if(layer == "f0")
                                        return line.MoveToX(dist);
                                    else
                                        return line.MoveToY(dist);
                                },
                                {"MergerData"})};
        // Fill histograms!
        ROOT::TThreadedObject<TH2D> hSP {ROOT::TNumSlots {node.GetNSlots()},
                                         "hSP",
                                         TString::Format("Side %.2f mm;X or Y [mm];Z [mm]", dist),
                                         xbins,
                                         xlims.first,
                                         xlims.second,
                                         zbins,
                                         zlims.first,
                                         zlims.second};
        std::map<int, ROOT::TThreadedObject<TH1D>> pxs, pzs;
        std::vector<int> idxs {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
        for(const auto& idx : idxs)
        {
            pxs.emplace(std::piecewise_construct, std::forward_as_tuple(idx),
                        std::forward_as_tuple(ROOT::TNumSlots {node.GetNSlots()}, TString::Format("px%d", idx),
                                              TString::Format("%.2f mm X or Y proj %d;X [mm]", dist, idx), xbins,
                                              xlims.first, xlims.second));
            pzs.emplace(std::piecewise_construct, std::forward_as_tuple(idx),
                        std::forward_as_tuple(ROOT::TNumSlots {node.GetNSlots()}, TString::Format("pz%d", idx),
                                              TString::Format("%.2f mm Z proj %d;Z [mm]", dist, idx), zbins,
                                              zlims.first, zlims.second));
        }
        // Create all slots 0 (blame your computer Iván)
        hSP.GetAtSlot(0)->GetEntries();
        for(auto m : {&pxs, &pzs})
            for(auto& [i, h] : *m)
                h.GetAtSlot(0)->GetEntries();

        node.ForeachSlot(
            [&](unsigned int slot, ActRoot::MergerData& data, ROOT::Math::XYZPointF& sp)
            {
                slot = (node.GetNSlots() - 1) - slot;
                if (layer == "f0")
                {
                    hSP.GetAtSlot(slot)->Fill(sp.Y(), sp.Z());
                }
                else
                {
                    hSP.GetAtSlot(slot)->Fill(sp.X(), sp.Z());
                }
                auto idx {data.fSilNs.front()};
                if(pxs.count(idx))
                {
                    if(layer == "f0")
                        pxs[idx].GetAtSlot(slot)->Fill(sp.Y());
                    else
                        pxs[idx].GetAtSlot(slot)->Fill(sp.X());
                    pzs[idx].GetAtSlot(slot)->Fill(sp.Z());
                }
            },
            {"MergerData", "NewSP"});

        // Write data
        f->cd();
        auto path {TString::Format("d_%.1f_mm/", dist)};
        auto* dir {f->mkdir(path)};
        dir->cd();
        hSP.Merge()->Write();
        // std::cout << "Entries in hSP: " << hSP.GetAtSlot(0)->GetEntries() << '\n';

        for(auto m : {&pxs, &pzs})
            for(auto& [i, h] : *m)
                h.Merge()->Write();

        f->cd();
    }
}
