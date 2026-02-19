#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActTypes.h"

#include "ROOT/RDF/HistoModels.hxx"
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


void checkKinematicLinesFront()
{
    PrettyStyle(true, true);
    std::string beam {"7Li"};
    std::string target {"d"};
    std::string light {"d"};

    // Get file output from pipe2
    TString filename {
        TString::Format("./Inputs/tree_ex_%s_%s_%s_filtered.root", beam.c_str(), target.c_str(), light.c_str())};
    ROOT::RDataFrame df {"Final_Tree", filename};

    // Ensure that are f0 events (but this trees were already filtered in Pipe1, so they should be)
    auto def = df.Filter([](ActRoot::MergerData& m)
                         { return m.fLight.GetNLayers() == 1 && m.fLight.GetLayer(0) == "f0"; }, {"MergerData"});

    // THETA LIGHT THETA HEAVY

    std::vector<std::string> listOfCuts11Libeam {"9Li", "11Li"};
    std::vector<std::string> listOfCuts7Libeam {"7Li", "8Li"};
    // Read cuts files for heavy particle, there is one per pad of f2 telescope (divide in 4 squares)
    ActRoot::CutsManager<std::string> cuts;
    // Read cuts on heavy particle
    for(const auto& recoil : listOfCuts11Libeam)
    {
        for(int s {0}; s < 4; s++)
        {
            auto cutfile {
                TString::Format("../../PostAnalysis/Cuts/pid_%s_f2_%d_11Li.root", recoil.c_str(), s)};
            cuts.ReadCut(TString::Format("%s_%d", recoil.c_str(), s).Data(), cutfile.Data());
        }
    }
    for(const auto& recoil : listOfCuts7Libeam)
    {
        for(int s {0}; s < 4; s++)
        {
            auto cutfile {
                TString::Format("../../PostAnalysis/Cuts/pid_%s_f2_%d_7Li.root", recoil.c_str(), s)};
            cuts.ReadCut(TString::Format("%s_%d", recoil.c_str(), s).Data(), cutfile.Data());
        }
    }

    // Filter de original dataframe with each cut checking whether the particle is inside it
    auto df7Li = def.Filter(
        [&cuts](ActRoot::MergerData& mer)
        {
            if(mer.fHeavy.GetNLayers() != 2)
                return false;
            else
            {
                auto n {std::to_string(mer.fHeavy.fNs[0])};
                auto key {TString::Format("7Li_%s", n.c_str())};
                return cuts.IsInside(key.Data(), mer.fHeavy.fEs[1], mer.fHeavy.fEs[0]);
            }
        },
        {"MergerData"});

    auto df8Li = def.Filter(
        [&cuts](ActRoot::MergerData& mer)
        {
            if(mer.fHeavy.GetNLayers() != 2)
                return false;
            else
            {
                auto n {std::to_string(mer.fHeavy.fNs[0])};
                auto key {TString::Format("8Li_%s", n.c_str())};
                return cuts.IsInside(key.Data(), mer.fHeavy.fEs[1], mer.fHeavy.fEs[0]);
            }
        },
        {"MergerData"});

    auto df9Li = def.Filter(
        [&cuts](ActRoot::MergerData& mer)
        {
            if(mer.fHeavy.GetNLayers() != 2)
                return false;
            else
            {
                auto n {std::to_string(mer.fHeavy.fNs[0])};
                auto key {TString::Format("9Li_%s", n.c_str())};
                return cuts.IsInside(key.Data(), mer.fHeavy.fEs[1], mer.fHeavy.fEs[0]);
            }
        },
        {"MergerData"});

    auto df11Li = def.Filter(
        [&cuts](ActRoot::MergerData& mer)
        {
            if(mer.fHeavy.GetNLayers() != 2)
                return false;
            else
            {
                auto n {std::to_string(mer.fHeavy.fNs[0])};
                auto key {TString::Format("11Li_%s", n.c_str())};
                return cuts.IsInside(key.Data(), mer.fHeavy.fEs[1], mer.fHeavy.fEs[0]);
            }
        },
        {"MergerData"});

    // Create histogram for ThetaLab - T3
    ROOT::RDF::TH2DModel Kin {
        "hKin", TString::Format("Kinematics for %s;#theta_{Lab} [#circ];E_{Vertex} [MeV]", light.c_str()),
        120,    0,
        60,     100,
        20,     45};
    ROOT::RDF::TH2DModel Kin_total {
        "hKin_total",
        TString::Format("Kinematics for %s ;#theta_{Lab} [#circ];E_{Vertex} [MeV]", light.c_str()),
        160,
        0,
        80,
        220,
        0,
        50};
    ROOT::RDF::TH1DModel Ex {"hEx",
                             TString::Format("Excitation energy for %s ;E_{x} [MeV];Counts / %.f keV", light.c_str(),
                                             (35. - (-10.)) / 300 * 1e3),
                             300, -10, 35};
    ROOT::RDF::TH2DModel ThetaLight_ThetaHeavy {
        "hThetaLight_ThetaHeavy",
        TString::Format("Theta light vs Theta heavy for %s;#theta_{Light} [#circ];#theta_{Heavy} [#circ]",
                        light.c_str()),
        120,
        0,
        60,
        120,
        0,
        30};

    auto hEx {def.Histo1D(Ex, "Ex")};
    auto hKin_total {def.Histo2D(Kin_total, "fThetaLight", "EVertex")};
    auto hThetaLight_ThetaHeavy {def.Histo2D(ThetaLight_ThetaHeavy, "fThetaLight", "fThetaHeavy")};

    auto hKin_total_7Li {df7Li.Histo2D(Kin_total, "fThetaLight", "EVertex")};
    hKin_total_7Li->SetTitle(
        TString::Format("Kinematics for %s with 7Li gate;#theta_{Lab} [#circ];E_{Vertex} [MeV]", light.c_str()));
    auto hEx_7Li {df7Li.Histo1D(Ex, "Ex")};
    hEx_7Li->SetTitle(TString::Format("Excitation energy for %s with 7Li gate;E_{x} [MeV];Counts / %.f keV",
                                      light.c_str(), (35. - (-10.)) / 300 * 1e3));
    auto hThetaLight_ThetaHeavy_7Li {df7Li.Histo2D(ThetaLight_ThetaHeavy, "fThetaLight", "fThetaHeavy")};
    hThetaLight_ThetaHeavy_7Li->SetTitle(
        TString::Format("Theta light vs Theta heavy for %s with 7Li gate;#theta_{Light} [#circ];#theta_{Heavy} [#circ]",
                        light.c_str()));
    auto hKin_total_8Li {df8Li.Histo2D(Kin_total, "fThetaLight", "EVertex")};
    hKin_total_8Li->SetTitle(
        TString::Format("Kinematics for %s with 8Li gate;#theta_{Lab} [#circ];E_{Vertex} [MeV]", light.c_str()));
    auto hEx_8Li {df8Li.Histo1D(Ex, "Ex")};
    hEx_8Li->SetTitle(TString::Format("Excitation energy for %s with 8Li gate;E_{x} [MeV];Counts / %.f keV",
                                      light.c_str(), (35. - (-10.)) / 300 * 1e3));
    auto hThetaLight_ThetaHeavy_8Li {df8Li.Histo2D(ThetaLight_ThetaHeavy, "fThetaLight", "fThetaHeavy")};
    hThetaLight_ThetaHeavy_8Li->SetTitle(
        TString::Format("Theta light vs Theta heavy for %s with 8Li gate;#theta_{Light} [#circ];#theta_{Heavy} [#circ]",
                        light.c_str()));
    auto hKin_total_9Li {df9Li.Histo2D(Kin_total, "fThetaLight", "EVertex")};
    hKin_total_9Li->SetTitle(
        TString::Format("Kinematics for %s with 9Li gate;#theta_{Lab} [#circ];E_{Vertex} [MeV]", light.c_str()));
    auto hEx_9Li {df9Li.Histo1D(Ex, "Ex")};
    hEx_9Li->SetTitle(TString::Format("Excitation energy for %s with 9Li gate;E_{x} [MeV];Counts / %.f keV",
                                      light.c_str(), (35. - (-10.)) / 300 * 1e3));
    auto hThetaLight_ThetaHeavy_9Li {df9Li.Histo2D(ThetaLight_ThetaHeavy, "fThetaLight", "fThetaHeavy")};
    hThetaLight_ThetaHeavy_9Li->SetTitle(
        TString::Format("Theta light vs Theta heavy for %s with 9Li gate;#theta_{Light} [#circ];#theta_{Heavy} [#circ]",
                        light.c_str()));
    auto hKin_total_11Li {df11Li.Histo2D(Kin_total, "fThetaLight", "EVertex")};
    hKin_total_11Li->SetTitle(
        TString::Format("Kinematics for %s with 11Li gate;#theta_{Lab} [#circ];E_{Vertex} [MeV]", light.c_str()));
    auto hEx_11Li {df11Li.Histo1D(Ex, "Ex")};
    hEx_11Li->SetTitle(TString::Format("Excitation energy for %s with 11Li gate;E_{x} [MeV];Counts / %.f keV",
                                       light.c_str(), (35. - (-10.)) / 300 * 1e3));
    auto hThetaLight_ThetaHeavy_11Li {df11Li.Histo2D(ThetaLight_ThetaHeavy, "fThetaLight", "fThetaHeavy")};
    hThetaLight_ThetaHeavy_11Li->SetTitle(TString::Format(
        "Theta light vs Theta heavy for %s with 11Li gate;#theta_{Light} [#circ];#theta_{Heavy} [#circ]",
        light.c_str()));

    // Theoretical kinematic lines
    ActPhysics::Particle pb {beam};
    ActPhysics::Particle pt {target};
    ActPhysics::Particle pl {light};
    double initialEnergy = 7.5; // Define initialEnergy appropriately
    ActPhysics::Kinematics kin_gs {pb, pt, pl, initialEnergy * pb.GetAMU(), 0};
    ActPhysics::Kinematics kin_4_1_ex {pb, pt, pl, initialEnergy * pb.GetAMU(), 4.1};

    auto* theo_gs {kin_gs.GetKinematicLine3()};
    auto* theo_4_1_ex {kin_4_1_ex.GetKinematicLine3()};
    auto* theo_Heavy_Light_gs {kin_gs.GetTheta3vs4Line()};

    auto* c0 {new TCanvas("c0", "Total plots")};
    c0->Divide(2, 2);
    c0->cd(1);
    hEx->DrawClone();
    c0->cd(2);
    hKin_total->DrawClone("colz");
    c0->cd(3);
    hThetaLight_ThetaHeavy->DrawClone("colz");
    theo_Heavy_Light_gs->Draw("same");

    auto* c3 {new TCanvas("c3", "Kinematics total with theo lines")};
    hKin_total->DrawClone("colz");
    // Draw kinematic lines
    theo_gs->Draw("same");
    theo_4_1_ex->Draw("same");

    auto* c4 {new TCanvas("c4", "Kinematic plots with heavy gates for 11Li beam")};
    c4->Divide(2, 2);
    c4->cd(1);
    hKin_total_9Li->DrawClone("colz");
    c4->cd(2);
    hKin_total_11Li->DrawClone("colz");
    c4->cd(3);
    hEx_9Li->DrawClone();
    c4->cd(4);
    hEx_11Li->DrawClone();

    auto* c5 {new TCanvas("c5", "Kinematic plots with heavy gates for 7Li beam")};
    c5->Divide(2, 2);
    c5->cd(1);
    hKin_total_7Li->DrawClone("colz");
    theo_gs->Draw("same");
    theo_4_1_ex->Draw("same");
    c5->cd(2);
    hKin_total_8Li->DrawClone("colz");
    c5->cd(3);
    hEx_7Li->DrawClone();
    c5->cd(4);
    hEx_8Li->DrawClone();

    auto* c6 {new TCanvas("c6", "Theta light vs Theta heavy with heavy gates")};
    c6->Divide(2, 2);
    c6->cd(1);
    hThetaLight_ThetaHeavy_9Li->DrawClone("colz");
    c6->cd(2);
    hThetaLight_ThetaHeavy_11Li->DrawClone("colz");
    c6->cd(3);
    hThetaLight_ThetaHeavy_7Li->DrawClone("colz");
    c6->cd(4);
    hThetaLight_ThetaHeavy_8Li->DrawClone("colz");
}