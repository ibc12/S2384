#ifndef Pipe2_ExM4_cxx
#define Pipe2_ExM4_cxx

#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActParticle.h"
#include "ActSRIM.h"
#include "ActSilData.h"
#include "ActSilMatrix.h"
#include "ActSilSpecs.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TH2.h"
#include "TH2D.h"
#include "TString.h"

#include <fstream>
#include <map>
#include <string>

#include "../../../PostAnalysis/HistConfig.h"
#include "../../../PrettyStyle.C"
#include "../Utils.h"

struct ExLevel
{
    double Eex;
    int color;
    std::string label;
};

void Pipe2_ExM4(const std::string& beam, const std::string& target, const std::string& light)
{
    // Get file from pipe1
    TString infile = TString::Format("./Outputs/PIDM4_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str());
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df("Final_Tree", infile.Data());

    std::string srimName {};
    if(light == "d")
        srimName = "2H";
    else if(light == "p")
        srimName = "1H";
    else if(light == "t")
        srimName = "3H";
    else if(light == "3He")
        srimName = "3He";
    else if(light == "4He")
        srimName = "4He";

    // Reconstruct ex with energy of protons or deuterons at vertex
    // For that we need the srim files n the gas for the light particles
    auto* srim {new ActPhysics::SRIM};
    srim->ReadTable(light, TString::Format("../../Calibrations/SRIM/%s_900mb_CF4_95-5.txt", srimName.c_str()).Data());
    srim->ReadTable(beam, TString::Format("../../Calibrations/SRIM/%s_900mb_CF4_95-5.txt", beam.c_str()).Data());
    srim->ReadTable("mylar", TString::Format("../../Calibrations/SRIM/%s_Mylar.txt", beam.c_str()).Data());

    ActPhysics::Particle pb {beam};
    ActPhysics::Particle pt {target};
    ActPhysics::Particle pl {light};

    // Initial energy
    double initialEnergy {7.558}; // meassured by operators; resolution of 0,19%
    initialEnergy = srim->Slow("mylar", initialEnergy * pb.GetAMU(), 0.0168);
    initialEnergy = srim->Slow(beam, initialEnergy, 60); // 60 mm of gas before the pad plane
    initialEnergy = initialEnergy / pb.GetAMU();         // back to amu units

    auto dfVertex =
        df.Define("EVertex",
                  [&](float TL, float eSil, std::string layer)
                  {
                      double ret {};
                      if(layer == "L1")
                          ret = srim->EvalEnergy(light, TL);
                      else if(layer == "l0" || layer == "r0" || layer == "f0")
                          ret = srim->EvalInitialEnergy(light, eSil, TL);
                      return ret;
                  },
                  {"TrackLength", "SilE", "Layer"})
            .Define("EBeam", [&](const ActRoot::TPCData& tpc)
                    { return srim->Slow(beam, initialEnergy * pb.GetAMU(), tpc.fRPs.front().X()); }, {"TPCData"});

    ActPhysics::Kinematics kin {pb, pt, pl, initialEnergy * pb.GetAMU()};
    std::vector<ActPhysics::Kinematics> vkins {dfVertex.GetNSlots()};
    for(auto& k : vkins)
        k = kin;
    dfVertex =
        dfVertex.DefineSlot("Ex",
                            [&](unsigned int slot, const double thetaLab, double EVertex, double EBeam)
                            {
                                vkins[slot].SetBeamEnergy(EBeam);
                                return vkins[slot].ReconstructExcitationEnergy(EVertex, (thetaLab)*TMath::DegToRad());
                            },
                            {"ThetaLab", "EVertex", "EBeam"});
    dfVertex = dfVertex.DefineSlot(
        "ThetaCM",
        [&](unsigned int slot, const double thetaLab, double EVertex, double EBeam)
        {
            vkins[slot].SetBeamEnergy(EBeam);
            return vkins[slot].ReconstructTheta3CMFromLab(EVertex, (thetaLab)*TMath::DegToRad()) * TMath::RadToDeg();
        },
        {"ThetaLab", "EVertex", "EBeam"});

    // Save dataframe in a .root file
    TString outfile = TString::Format("./Outputs/ExM4_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str());
    dfVertex.Snapshot("Final_Tree", outfile);
    std::cout << "Saving Final_Tree in " << outfile << " from Pipe2_ExM4" << '\n';

    // ---- Separamos por layer: L1 vs el resto (l0, r0, f0) ----
    auto dfL1 = dfVertex.Filter([](const std::string& layer) { return layer == "L1"; }, {"Layer"});
    auto dfOthers = dfVertex.Filter([](const std::string& layer) { return layer != "L1"; }, {"Layer"});

    // Histogramas de Ex para cada subconjunto
    auto hExL1 {dfL1.Histo1D(HistConfig::Ex, "Ex")};
    auto hExOthers {dfOthers.Histo1D(HistConfig::Ex, "Ex")};

    // Histogramas de cinematica para cada subconjunto
    auto hkinL1 {dfL1.Histo2D(HistConfig::Kin, "ThetaLab", "EVertex")};
    auto hkinOthers {dfOthers.Histo2D(HistConfig::Kin, "ThetaLab", "EVertex")};

    // --- Plots de Ex ---
    auto c1 = new TCanvas("cExM4_L1", "Excitation energy Multiplicity 4 - L1", 800, 600);
    c1->cd();
    hExL1->SetTitle("Ex, layer L1");
    hExL1->DrawClone();

    auto c2 = new TCanvas("cExM4_Others", "Excitation energy Multiplicity 4 - Other layers", 800, 600);
    c2->cd();
    hExOthers->SetTitle("Ex, layers l0/r0/f0");
    hExOthers->DrawClone();

    // --- Plots de cinematica ---
    auto c3 = new TCanvas("cKinM4_L1", "Kinematics Multiplicity 4 - L1", 800, 600);
    c3->cd();
    hkinL1->SetTitle("Kinematics, layer L1");
    hkinL1->DrawClone("colz");

    auto c4 = new TCanvas("cKinM4_Others", "Kinematics Multiplicity 4 - Other layers", 800, 600);
    c4->cd();
    hkinOthers->SetTitle("Kinematics, layers l0/r0/f0");
    hkinOthers->DrawClone("colz");

    // theo kin: una linea por cada nivel de excitacion, con su color y su entrada en la leyenda
    // Se dibuja sobre ambos canvases de cinematica (L1 y otras layers)

    std::vector<ExLevel> levels {
        {0.0, kBlack, "g.s."},
        {4.63, kRed, "4.63 MeV"},
        {6.53, kBlue, "6.53 MeV"},
        {7.1, kGreen + 2, "7.1 MeV"},
    };

    auto* legendL1 = new TLegend(0.65, 0.65, 0.88, 0.88);
    legendL1->SetBorderSize(0);
    legendL1->SetFillStyle(0);

    auto* legendOthers = new TLegend(0.65, 0.65, 0.88, 0.88);
    legendOthers->SetBorderSize(0);
    legendOthers->SetFillStyle(0);

    // Guardamos los objetos Kinematics para que no mueran al salir del bucle
    static std::vector<ActPhysics::Kinematics> kinLevels;
    kinLevels.clear();
    kinLevels.reserve(levels.size());

    for(const auto& level : levels)
    {
        kinLevels.emplace_back(pb, pt, pl, initialEnergy * pb.GetAMU(), level.Eex);
        auto* theo {kinLevels.back().GetKinematicLine3()};
        theo->SetLineColor(level.color);
        theo->SetLineWidth(2);

        c3->cd();
        auto* theoL1 {(TGraph*)theo->DrawClone("same")};
        legendL1->AddEntry(theoL1, level.label.c_str(), "l");

        c4->cd();
        auto* theoOthers {(TGraph*)theo->DrawClone("same")};
        legendOthers->AddEntry(theoOthers, level.label.c_str(), "l");
    }

    c3->cd();
    legendL1->Draw();

    c4->cd();
    legendOthers->Draw();

    std::ofstream out(TString::Format("./Outputs/Pipe2_ExM4_%s_ex.dat", light.c_str()).Data());
    dfL1.Foreach(
        [&](ActRoot::MergerData& m, double ex)
        {
            if(ex > 2)
                m.Stream(out);
        },
        {"MergerData", "Ex"});
    out.close();
}
#endif