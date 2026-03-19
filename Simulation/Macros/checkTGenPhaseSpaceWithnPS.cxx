#include "ActColors.h"
#include "ActCrossSection.h"
#include "ActCutsManager.h"
#include "ActDecayGenerator.h"
#include "ActKinematicGenerator.h"
#include "ActKinematics.h"
#include "ActLine.h"
#include "ActParticle.h"
#include "ActRunner.h"
#include "ActSRIM.h"
#include "ActSilData.h"
#include "ActSilSpecs.h"
#include "ActTPCParameters.h"
#include "ActUtils.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TEfficiency.h"
#include "TF1.h"
#include "TFile.h"
#include "TGenPhaseSpace.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TPolyLine3D.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>


void checkTGenPhaseSpaceWithnPS()
{
    std::string beam {"11Li"};
    std::string target {"d"};
    std::string light {"d"};
    std::string heavy {"11Li"};
    int neutronPS {1};       // number of neutrons in final state
    int protonPS {0};        // number of protons in final state
    double Tbeam {11 * 7.5}; // MeV
    double Ex {0};           // MeV

    auto* kin {new ActPhysics::Kinematics {beam, target, light, heavy, Tbeam, Ex}};
    auto* kinGen {new ActSim::KinematicGenerator {beam, target, light, heavy, protonPS, neutronPS}};

    // Sample kinematics generator
    kinGen->SetBeamAndExEnergies(Tbeam, 0);
    double weight = kinGen->Generate();
    if(neutronPS == 0 && protonPS == 0)
    {
        weight = 1;
    }
    kin = kinGen->GetBinaryKinematics();
    // Get Lorenzt vector of products
    auto LorenztVector3 = kinGen->GetLorentzVector(0);
    auto LorenztVector4 = kinGen->GetLorentzVector(1);
    // Get angles
    double theta3Lab = LorenztVector3->Theta();
    double phi3Lab = LorenztVector3->Phi();
    double T3Lab = LorenztVector3->E() - LorenztVector3->M();
    // Save without resolution
    double theta3LabSampled = theta3Lab;

    double theta3CMBefore = kin->ReconstructTheta3CMFromLab(T3Lab, theta3LabSampled) *
                            TMath::RadToDeg(); // this is in deg, because of xs sampling in other case
    double theta3CM = kin->ReconstructTheta3CMFromLab(T3Lab, theta3Lab);
    double phi3CM = phi3Lab;

    // Heavy
    double theta4Lab = LorenztVector4->Theta();
    double phi4Lab = LorenztVector4->Phi();
    double T4Lab = LorenztVector4->E() - LorenztVector4->M();

    std::cout << "Mass of particle 4: " << LorenztVector4->M() / 931.5 << std::endl;


    // Get root file from Outputs to plot Ex from PS
    auto filename {"../Outputs/11Li/2H_2H_TRIUMF_Eex_0.000_nPS_2_pPS_0.root"};
    // auto filename {TString::Format("./Outputs/tree_pid_11Li_d_p.root")};
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"SimulationTTree", filename};

    auto hEx = df.Histo1D({"hEx", "Excitation energy from PS;E_{x} [MeV];Counts", 100, -1, 10}, "EexGateHeavy");
    auto* c = new TCanvas();
    hEx->DrawClone();
}