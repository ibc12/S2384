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
#include "TLegend.h"
#include "TString.h"

#include <map>
#include <string>

#include "../../PostAnalysis/HistConfig.h"
#include "./Utils.h"

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
                  [&](float TL, float eSil)
                  {
                      double ret = srim->EvalInitialEnergy(light, eSil, TL);
                      // else // L1 trigger
                      //     ret = srim->EvalEnergy("p", d.fLight.fTL);
                      return ret;
                  },
                  {"TrackLength", "SilESelectedParticle"})
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

    // Plot Ex and kin
    auto hEx {dfVertex.Histo1D(HistConfig::Ex, "Ex")};
    auto hkin {dfVertex.Histo2D(HistConfig::Kin, "ThetaLab", "EVertex")};

    auto c1 = new TCanvas("cExM4_0", "Excitation energy Multiplicity 4", 800, 600);
    c1->cd();
    hEx->DrawClone();
    auto c2 = new TCanvas("cExM4_1", "Kinematics Multiplicity 4", 800, 600);
    c2->cd();
    hkin->DrawClone("colz");
    // theo kin: una linea por cada nivel de excitacion, con su color y su entrada en la leyenda

    std::vector<ExLevel> levels {
        {0.0, kBlack, "g.s."},
        {4.63, kRed, "4.63 MeV"},
        {6.53, kBlue, "6.53 MeV"},
        {7.1, kGreen + 2, "7.1 MeV"},
    };

    auto* legend = new TLegend(0.65, 0.65, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);

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
        theo->Draw("same");
        legend->AddEntry(theo, level.label.c_str(), "l");
    }

    legend->Draw();
}
#endif