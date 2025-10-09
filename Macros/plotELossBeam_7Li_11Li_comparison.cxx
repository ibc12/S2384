#include "ActDataManager.h"
#include "ActTPCData.h"
#include "ActSilData.h"
#include "ActModularData.h"
#include "ActMergerData.h"
#include <ROOT/RDataFrame.hxx>
#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"
#include "TLegend.h"

void plotELossBeam_7Li_11Li_comparison()
{
    // get 1 run from 7Li and other from 11Li
    ActRoot::DataManager dataman7Li{"../configs/data_7Li.conf", ActRoot::ModeType::EReadSilMod};
    auto chain7Li{dataman7Li.GetChain()};
    auto friend1_7Li{dataman7Li.GetChain(ActRoot::ModeType::EMerge)};
    chain7Li->AddFriend(friend1_7Li.get());
    auto friend2_7Li{dataman7Li.GetChain(ActRoot::ModeType::EFilter)};
    chain7Li->AddFriend(friend2_7Li.get());
    ROOT::RDataFrame df7Li{*chain7Li};

    ActRoot::DataManager dataman11Li{"../configs/data.conf", ActRoot::ModeType::EReadSilMod};
    auto chain11Li{dataman11Li.GetChain()};
    auto friend1_11Li{dataman11Li.GetChain(ActRoot::ModeType::EMerge)};
    chain11Li->AddFriend(friend1_11Li.get());
    auto friend2_11Li{dataman11Li.GetChain(ActRoot::ModeType::EFilter)};
    chain11Li->AddFriend(friend2_11Li.get());
    ROOT::RDataFrame df11Li{*chain11Li};

    // Create columns for E_Loss and E_Beam
    auto df7Li_gated = df7Li.Filter([](ActRoot::MergerData &m)
                                    { return m.fBeamIdx != -1; }, {"MergerData"});
    auto df11Li_gated = df11Li.Filter([](ActRoot::MergerData &m)
                                      { return m.fBeamIdx != -1; }, {"MergerData"});

    auto df7Li_final = df7Li_gated.Define("E_Loss", [](ActRoot::MergerData &m, ActRoot::TPCData &d)
                                          {int beamIdx {m.fBeamIdx}; 
                                             auto voxelsBeam {d.fClusters[beamIdx].GetRefToVoxels()};
                                             double DE {0};
                                             for (auto v : voxelsBeam)
                                             {
                                                if(v.GetPosition().X() < 10)
                                                {
                                                    DE += v.GetCharge();
                                                }
                                             }
                                             return DE; }, {"MergerData", "TPCData"})
                           .Define("E_Beam", [](ActRoot::MergerData &m, ActRoot::TPCData &d)
                                   { int beamIdx {m.fBeamIdx}; 
                                             auto voxelsBeam {d.fClusters[beamIdx].GetRefToVoxels()};
                                             double E {0};
                                             for (auto v : voxelsBeam)
                                             {
                                                E += v.GetCharge();
                                             }
                                             return E; }, {"MergerData", "TPCData"});

    auto df11Li_final = df11Li_gated.Define("E_Loss", [](ActRoot::MergerData &m, ActRoot::TPCData &d)
                                            { int beamIdx {m.fBeamIdx}; 
                                             auto voxelsBeam {d.fClusters[beamIdx].GetRefToVoxels()};
                                             double DE {};
                                             for (auto v : voxelsBeam)
                                             {
                                                if(v.GetPosition().X() < 10)
                                                {
                                                    DE += v.GetCharge();
                                                }
                                             }
                                             return DE; }, {"MergerData", "TPCData"})
                            .Define("E_Beam", [](ActRoot::MergerData &m, ActRoot::TPCData &d)
                                    { int beamIdx {m.fBeamIdx}; 
                                             auto voxelsBeam {d.fClusters[beamIdx].GetRefToVoxels()};
                                             double E {};
                                             for (auto v : voxelsBeam)
                                             {
                                                E += v.GetCharge();
                                             }
                                             return E; }, {"MergerData", "TPCData"});

    auto hDE_E_7Li{df7Li_final.Histo2D({"hDE_E_7Li", "DE vs E;DE [u.a.];E [u.a.]", 100, 0, 1e5, 100, 0, 1e5}, "E_Beam", "E_Loss")};
    auto hDE_E_11Li{df11Li_final.Histo2D({"hDE_E_11Li", "DE vs E;DE [u.a.];E [u.a.]", 100, 0, 1e5, 100, 0, 1e5}, "E_Beam", "E_Loss")};

    // Plot
    auto c{new TCanvas("c")};
    hDE_E_7Li->DrawClone("colz");
    hDE_E_11Li->DrawClone("same");

    // Extract vectors from RDataFrame
    std::vector<double> x7, y7, x11, y11;

    // For 7Li
    df7Li_final.Foreach([&](double Ebeam, double Eloss){ 
        x7.push_back(Ebeam); 
        y7.push_back(Eloss); 
    }, {"E_Beam","E_Loss"});

    // For 11Li
    df11Li_final.Foreach([&](double Ebeam, double Eloss){ 
        x11.push_back(Ebeam); 
        y11.push_back(Eloss); 
    }, {"E_Beam","E_Loss"});

    // Create TGraphs
    auto g7Li = new TGraph(x7.size(), x7.data(), y7.data());
    g7Li->SetMarkerStyle(20);
    g7Li->SetMarkerColor(kRed);
    g7Li->SetTitle("E_{Loss} vs E_{Beam};E_{Beam} [a.u.];E_{Loss} [a.u.]");

    auto g11Li = new TGraph(x11.size(), x11.data(), y11.data());
    g11Li->SetMarkerStyle(21);
    g11Li->SetMarkerColor(kBlue);

    // Draw
    auto c1 = new TCanvas("c1", "E_Loss vs E_Beam", 800,600);
    g7Li->Draw("AP");
    g11Li->Draw("P same");

    // Legend
    auto leg = new TLegend(0.65,0.75,0.88,0.88);
    leg->AddEntry(g7Li, "^{7}Li","p");
    leg->AddEntry(g11Li, "^{11}Li","p");
    leg->Draw();
}