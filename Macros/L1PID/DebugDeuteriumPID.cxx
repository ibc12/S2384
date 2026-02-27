#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TH2.h"
#include "TH2D.h"
#include "TString.h"

#include <fstream>
#include <map>
#include <string>

#include "../../PrettyStyle.C"

void DebugDeuteriumPID()
{
    PrettyStyle(true, true);
    // Get all data for 7Li (there is no triton there, so easier to see deuterium)
    std::string beam {"7Li"};
    std::string target {"d"};
    std::string light {"d"};

    std::string dataconf {};
    if(beam == "11Li")
        dataconf = "../../configs/data_11Li.conf";
    else if(beam == "7Li")
        dataconf = "../../configs/data_7Li.conf";
    else
        throw std::runtime_error("Beam cannot differ from 11Li or 7Li");

    // Read data
    ActRoot::DataManager dataman {dataconf, ActRoot::ModeType::EMerge};
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(chain2.get());

    // RDataFrame
    // ROOT::EnableImplicitMT();
    ROOT::RDataFrame dforigin {*chain};

    // Filter GATCONF == L1
    auto df = dforigin.Filter([](ActRoot::ModularData& mod) { return mod.Get("GATCONF") == 8; }, {"ModularData"});
    // Get the cut for PID
    ActRoot::CutsManager<std::string> cuts;
    cuts.ReadCut("l1",
                 TString::Format("../../PostAnalysis/Cuts/pid_%s_l1_%s.root", light.c_str(), beam.c_str()).Data());
    cuts.ReadCut("l1_theta",
                 TString::Format("./Cuts/%s_gs_ThetaVSq_%s.root", light.c_str(), beam.c_str()).Data());

    auto dfL1Cut = df.Filter(
        [&cuts](const ActRoot::MergerData& m)
        {
            if(cuts.GetCut("l1"))
                return cuts.IsInside("l1", m.fLight.fRawTL, m.fLight.fQtotal);
            else
                return false;
        },
        {"MergerData"});

    // Fill histogram TL vs Q and ThetaLight vs Q for both dfs
    auto hL1_track = new TH2D("hL1_track", "L1 TL vs Q;TL [mm];Q_{total}", 240, 0, 120, 2000, 0, 300000);
    auto hL1_track_cut = new TH2D("hL1_track_cut", "L1 TL vs Q (cut);TL [mm];Q_{total}", 240, 0, 120, 2000, 0, 300000);

    auto hL1_theta =
        new TH2D("hL1_theta", "L1 #theta_{light} vs Q;#theta_{light} [deg];Q_{total}", 270, 0, 180, 2000, 0, 300000);
    auto hL1_theta_cut = new TH2D("hL1_theta_cut", "L1 #theta_{light} vs Q (cut);#theta_{light} [deg];Q_{total}", 270,
                                  0, 180, 2000, 0, 300000);
    df.Foreach(
        [&hL1_track, &hL1_theta](const ActRoot::MergerData& m)
        {
            hL1_track->Fill(m.fLight.fRawTL, m.fLight.fQtotal);
            hL1_theta->Fill(m.fThetaLight, m.fLight.fQtotal);
        },
        {"MergerData"});

    dfL1Cut.Foreach(
        [&hL1_track_cut, &hL1_theta_cut](const ActRoot::MergerData& m)
        {
            hL1_track_cut->Fill(m.fLight.fRawTL, m.fLight.fQtotal);
            hL1_theta_cut->Fill(m.fThetaLight, m.fLight.fQtotal);
        },
        {"MergerData"});

    // Save events inside cut Theta vs Q
    std::ofstream outFile(TString::Format("./Outputs/events_inside_cut_%s_l1Theta_%s.dat", light.c_str(), beam.c_str()).Data());
    dfL1Cut.Foreach(
        [&outFile, &cuts](const ActRoot::MergerData& m)        {
            if(cuts.GetCut("l1_theta"))
                if(cuts.IsInside("l1_theta", m.fThetaLight, m.fLight.fQtotal))
                m.Stream(outFile);
        },
        {"MergerData"});
        outFile.close();

    // Draw the histograms
    auto* c = new TCanvas("c", "Deuterium PID", 1200, 600);
    c->Divide(2, 2);
    c->cd(1);
    hL1_track->Draw("COLZ");
    cuts.GetCut("l1")->Draw("same");
    c->cd(2);
    hL1_track_cut->Draw("COLZ");
    c->cd(3);
    hL1_theta->Draw("COLZ");
    c->cd(4);
    hL1_theta_cut->Draw("COLZ");
    cuts.GetCut("l1_theta")->Draw("same");
}