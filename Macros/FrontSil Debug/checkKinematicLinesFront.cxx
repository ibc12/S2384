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
    std::string beam {"11Li"};
    std::string target {"d"};
    std::string light {"t"};

    // Get file output from pipe2
    TString filename {
        TString::Format("./Inputs/tree_ex_%s_%s_%s_filtered.root", beam.c_str(), target.c_str(), light.c_str())};
    ROOT::RDataFrame df {"Final_Tree", filename};

    // Ensure that are f0 events (but this trees were already filtered in Pipe1, so they should be)
    auto def = df.Filter([](ActRoot::MergerData& m)
                         { return m.fLight.GetNLayers() == 1 && m.fLight.GetLayer(0) == "f0"; }, {"MergerData"});

    auto defFiltered_gs = def.Filter(
        [](double& ex)
        {
            if(ex < 1 && ex > -1)
                return true;
            else
                return false;
        },
        {"Ex"});

    auto defFiltered_4_1_ex = def.Filter(
        [](double& ex)
        {
            if(ex < 4.4 && ex > 3.9)
                return true;
            else
                return false;
        },
        {"Ex"});

    // Create histogram for ThetaLab - T3
    ROOT::RDF::TH2DModel Kin {"hKin", TString::Format("Kinematics for %s;#theta_{Lab} [#circ];E_{Vertex} [MeV]", light.c_str()), 120, 0, 60, 100, 20, 45};
    ROOT::RDF::TH2DModel Kin_total {
        "hKin_total", TString::Format("Kinematics for %s ;#theta_{Lab} [#circ];E_{Vertex} [MeV]", light.c_str()), 160, 0, 80, 220, 0, 45};
    ROOT::RDF::TH1DModel Ex {
        "hEx", TString::Format("Excitation energy for %s ;E_{x} [MeV];Counts / %.f keV", light.c_str(), (35. - (-10.)) / 300 * 1e3), 300, -10,
        35};

    auto hEx {def.Histo1D(Ex, "Ex")};

    auto hKin_gs {defFiltered_gs.Histo2D(Kin, "fThetaLight", "EVertex")};
    auto hKin_4_1_ex {defFiltered_4_1_ex.Histo2D(Kin, "fThetaLight", "EVertex")};
    auto hKin_total {def.Histo2D(Kin_total, "fThetaLight", "EVertex")};

    // Theoretical kinematic lines
    ActPhysics::Particle pb {beam};
    ActPhysics::Particle pt {target};
    ActPhysics::Particle pl {light};
    double initialEnergy = 7.5; // Define initialEnergy appropriately
    ActPhysics::Kinematics kin_gs {pb, pt, pl, initialEnergy * pb.GetAMU(), 0};
    ActPhysics::Kinematics kin_4_1_ex {pb, pt, pl, initialEnergy * pb.GetAMU(), 4.1};

    auto* theo_gs {kin_gs.GetKinematicLine3()};
    auto* theo_4_1_ex {kin_4_1_ex.GetKinematicLine3()};

    auto* c0 {new TCanvas("c0", "Excitation energy")};
    hEx->DrawClone();

    auto* c1 {new TCanvas("c1", "Kinematics zoom")};
    hKin_gs->DrawClone("colz");
    hKin_4_1_ex->DrawClone("same");
    theo_gs->Draw("same");
    theo_4_1_ex->Draw("same");

    auto* c2 {new TCanvas("c2", "Kinematics total")};
    hKin_total->DrawClone("colz");

    auto* c3 {new TCanvas("c3", "Kinematics total with theo lines")};
    hKin_total->DrawClone("colz");
    // Draw kinematic lines
    theo_gs->Draw("same");
    theo_4_1_ex->Draw("same");
}