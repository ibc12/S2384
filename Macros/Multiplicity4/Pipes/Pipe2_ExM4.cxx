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

#include <map>
#include <string>

#include "../../../PrettyStyle.C"
#include "../../../PostAnalysis/HistConfig.h"
#include "../Utils.h"

void Pipe2_ExM4(const std::string& beam, const std::string& target, const std::string& light)
{
    // Get file from pipe1
    TString infile = TString::Format("./Outputs/PIDM4_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str());
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df("Final_Tree", infile.Data());

    // Reconstruct ex with energy of protons or deuterons at vertex
    // For that we need the srim files n the gas for the light particles
    auto* srim {new ActPhysics::SRIM};
    srim->ReadTable("p", TString::Format("../../Calibrations/SRIM/1H_900mb_CF4_95-5.txt").Data());
    srim->ReadTable("d", TString::Format("../../Calibrations/SRIM/2H_900mb_CF4_95-5.txt").Data());
    srim->ReadTable("7Li", TString::Format("../../Calibrations/SRIM/7Li_900mb_CF4_95-5.txt").Data());
    srim->ReadTable("mylar", TString::Format("../../Calibrations/SRIM/7Li_Mylar.txt").Data());

    ActPhysics::Particle pb {"7Li"};
    ActPhysics::Particle pt {"d"};
    ActPhysics::Particle pl {"p"};

    // Initial energy
    double initialEnergy {7.558}; // meassured by operators; resolution of 0,19%
    initialEnergy = srim->Slow("mylar", initialEnergy * pb.GetAMU(), 0.0168);
    initialEnergy = srim->Slow("7Li", initialEnergy, 60); // 60 mm of gas before the pad plane
    initialEnergy = initialEnergy / pb.GetAMU();          // back to amu units

    auto dfVertex =
        df.Define("EVertex",
                  [&](float distGas, float eSil)
                  {
                      double ret = srim->EvalInitialEnergy("p", eSil, distGas);
                      // else // L1 trigger
                      //     ret = srim->EvalEnergy("p", d.fLight.fTL);
                      return ret;
                  },
                  {"DistanceInGas", "SilESelectedParticle"})
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
}
#endif