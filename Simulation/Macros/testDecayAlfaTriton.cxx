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

void testDecayAlfaTriton()
{
    // Create a 12Li particle with a specific momentum
    std::string beam = "7Li";
    std::string target = "2H";
    std::string light = "1H";
    std::string heavy = "8Li";
    int pPS = 0;
    int nPS = 1;
    std::string decay1 = "4He";
    std::string decay2 = "3H";
    double TBeam = 7 * 7.558; // MeV
    double Ex = 4.63;


    auto* kin {new ActPhysics::Kinematics {beam, target, light, heavy, TBeam, Ex}};
    auto* kinGen {new ActSim::KinematicGenerator {beam, target, light, heavy, pPS, nPS}};

    auto hThetaProtonLab = new TH1D("hThetaProtonLab", "Proton Lab Angle", 100, 0, 180);
    auto hThetaHeavyLab = new TH1D("hThetaHeavyLab", "Heavy Lab Angle", 100, 0, 180);
    auto hThetaAlfaLab = new TH1D("hThetaAlfaLab", "Alfa Lab Angle", 100, 0, 180);
    auto hThetaTritonLab = new TH1D("hThetaTritonLab", "Triton Lab Angle", 100, 0, 180);

    auto hEnergyProtonLab = new TH1D("hEnergyProtonLab", "Proton Lab Energy", 100, 0, 20);
    auto hEnergyHeavyLab = new TH1D("hEnergyHeavyLab", "Heavy Lab Energy", 100, 0, 50);
    auto hEnergyAlfaLab = new TH1D("hEnergyAlfaLab", "Alfa Lab Energy", 100, 0, 50);
    auto hEnergyTritonLab = new TH1D("hEnergyTritonLab", "Triton Lab Energy", 100, 0, 50);

    auto hEnergyDifference = new TH1D("hEnergyDifference", "Energy Difference", 100, -5, 5);

    for(int i = 0; i < 100000; ++i)
    {
        double T3Lab {};
        double T4Lab {};
        double theta3Lab {};
        double theta3LabSampled {};
        double phi3Lab {};
        double theta4Lab {};
        double phi4Lab {};
        double theta3CM {};
        double phi3CM {};
        double theta3CMBefore {-1};
        double weight {1.};
        double EnergyDifference = 0.0;

        // Sample kinematics generator
        kinGen->SetBeamAndExEnergies(TBeam, Ex);
        weight = kinGen->Generate();
        if(nPS == 0 && pPS == 0)
        {
            weight = 1;
        }
        kin = kinGen->GetBinaryKinematics();
        // Get Lorenzt vector of products
        auto LorenztVector3 = kinGen->GetLorentzVector(0);
        auto LorenztVector4 = kinGen->GetLorentzVector(1);
        // Get angles
        theta3Lab = LorenztVector3->Theta();
        phi3Lab = LorenztVector3->Phi();
        T3Lab = LorenztVector3->E() - LorenztVector3->M();

        theta3CMBefore = kin->ReconstructTheta3CMFromLab(T3Lab, theta3LabSampled) *
                         TMath::RadToDeg(); // this is in deg, because of xs sampling in other case
        theta3CM = kin->ReconstructTheta3CMFromLab(T3Lab, theta3Lab);
        phi3CM = phi3Lab;

        // Heavy
        theta4Lab = LorenztVector4->Theta();
        phi4Lab = LorenztVector4->Phi();
        T4Lab = LorenztVector4->E() - LorenztVector4->M();

        // Decay heavy to alfa + triton
        // 1. Get particles to decay
        auto heavyParticle = ActPhysics::Particle("7Li");
        heavyParticle.SetEx(Ex);
        auto decay1Particle = ActPhysics::Particle(decay1);
        auto decay2Particle = ActPhysics::Particle(decay2);
        auto decayGen7Li = new ActSim::DecayGenerator(heavyParticle, decay1Particle, decay2Particle);
        decayGen7Li->SetDecay(T4Lab, theta4Lab, phi4Lab);
        decayGen7Li->Generate();
        auto decayLorenztVector_alfa = decayGen7Li->GetLorentzVector(0);
        auto decayLorenztVector_triton = decayGen7Li->GetLorentzVector(1);

        auto decayTheta1Lab = decayLorenztVector_alfa->Theta();
        auto decayTLab_alfa = decayLorenztVector_alfa->E() - decayLorenztVector_alfa->M();

        auto decayTheta2Lab = decayLorenztVector_triton->Theta();
        auto decayTLab_triton = decayLorenztVector_triton->E() - decayLorenztVector_triton->M();

        EnergyDifference = decayTLab_alfa + decayTLab_triton - T4Lab;

        // Fill histograms
        hThetaProtonLab->Fill(theta3Lab * TMath::RadToDeg(), weight);
        hThetaHeavyLab->Fill(theta4Lab * TMath::RadToDeg(), weight);
        hThetaAlfaLab->Fill(decayTheta1Lab * TMath::RadToDeg(), weight);
        hThetaTritonLab->Fill(decayTheta2Lab * TMath::RadToDeg(), weight);
        hEnergyProtonLab->Fill(T3Lab, weight);
        hEnergyHeavyLab->Fill(T4Lab, weight);
        hEnergyAlfaLab->Fill(decayTLab_alfa, weight);
        hEnergyTritonLab->Fill(decayTLab_triton, weight);
        hEnergyDifference->Fill(EnergyDifference, weight);
    }

    auto* c = new TCanvas("c", "c", 1200, 800);
    c->Divide(2, 2);
    c->cd(1);
    hThetaProtonLab->DrawClone("hist");
    c->cd(2);
    hThetaHeavyLab->DrawClone("hist");
    c->cd(3);
    hThetaAlfaLab->DrawClone("hist");
    c->cd(4);
    hThetaTritonLab->DrawClone("hist");

    auto* c2 = new TCanvas("c2", "c2", 1200, 800);
    c2->Divide(2, 2);
    c2->cd(1);
    hEnergyProtonLab->DrawClone("hist");
    c2->cd(2);
    hEnergyHeavyLab->DrawClone("hist");
    c2->cd(3);
    hEnergyAlfaLab->DrawClone("hist");
    c2->cd(4);
    hEnergyTritonLab->DrawClone("hist");

    auto* c3 = new TCanvas("c3", "c3", 800, 600);
    hEnergyDifference->DrawClone("hist");
}