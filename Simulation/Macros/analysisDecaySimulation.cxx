
#include "ActDecayGenerator.h"
#include "ActKinematicGenerator.h"
#include "ActKinematics.h"
#include "ActParticle.h"
#include "ActSRIM.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TH1.h"
#include "TMath.h"
#include "TROOT.h"

#include <cmath>
#include <iostream>
#include <string>


void analysisDecaySimulation()
{
    ROOT::RDataFrame df {"SimulationTTree",
                         "../Outputs/7Li/2H_1H_TRIUMF_Eex_7.100_nPS_1_pPS_0_decay_democratic4Body.root"};

    auto dfTriton = df.Filter("tritonLayerCode != 0", "Select events where triton is detected");
    auto dfAlfa = df.Filter("alfaLayerCode != 0", "Select events where alfa is detected");
    auto dfAlfaTriton = df.Filter("gateAlfaTriton != 0", "Select events where both alfa and triton are detected");
    auto dfOnlyTriton =
        df.Filter("tritonLayerCode != 0 && gateAlfaTriton == 0", "Select events where triton is detected but not alfa");
    auto dfOnlyAlfa =
        df.Filter("alfaLayerCode != 0 && gateAlfaTriton == 0", "Select events where alfa is detected but not triton");

    // Start by reading some SRIM tables to get the ranges of the particles in the gas
    auto* srim {new ActPhysics::SRIM};
    srim->ReadTable("alfaGas", "../../Calibrations/SRIM/4He_900mb_CF4_95-5.txt");
    srim->ReadTable("tritonGas", "../../Calibrations/SRIM/3H_900mb_CF4_95-5.txt");

    auto minRangeAlfa {srim->EvalRange("alfaGas", 6.0)};     // MeV
    auto minRangeTriton {srim->EvalRange("tritonGas", 4.0)}; // MeV

    std::cout << "Minimum range for alfa: " << minRangeAlfa << " mm" << std::endl;
    std::cout << "Minimum range for triton: " << minRangeTriton << " mm" << std::endl;

    // Let's do the kinematic plot for all the cases
    auto hKin_protonAll =
        df.Histo2D({"hKin_protonAll", "Proton Kinematics;Theta [deg];Energy [MeV]", 100, 0, 180, 100, 0, 20},
                   "theta3Lab", "EVertex");
    auto hKin_protonOnlyTriton = dfOnlyTriton.Histo2D(
        {"hKin_protonOnlyTriton", "Proton Kinematics (only triton);Theta [deg];Energy [MeV]", 100, 0, 180, 100, 0, 20},
        "theta3Lab", "EVertex");
    auto hKin_protonOnlyAlfa = dfOnlyAlfa.Histo2D(
        {"hKin_protonOnlyAlfa", "Proton Kinematics (only alfa);Theta [deg];Energy [MeV]", 100, 0, 180, 100, 0, 20},
        "theta3Lab", "EVertex");
    auto hKin_protonAlfaTriton = dfAlfaTriton.Histo2D(
        {"hKin_protonAlfaTriton", "Proton Kinematics (alfa+triton);Theta [deg];Energy [MeV]", 100, 0, 180, 100, 0, 20},
        "theta3Lab", "EVertex");

    // Now kinematic plot for triton and alfa
    auto hKin_tritonAll =
        df.Histo2D({"hKin_tritonAll", "Triton Kinematics;Theta [deg];Energy [MeV]", 100, 0, 180, 100, 0, 40},
                   "tritonAngle", "tritonEnergy");
    auto hKin_alfaAll =
        df.Histo2D({"hKin_alfaAll", "Alfa Kinematics;Theta [deg];Energy [MeV]", 100, 0, 180, 100, 0, 40}, "alfaAngle",
                   "alfaEnergy");

    auto hKin_tritonAlfaTriton = dfAlfaTriton.Histo2D(
        {"hKin_tritonAlfaTriton", "Triton Kinematics (alfa+triton);Theta [deg];Energy [MeV]", 100, 0, 180, 100, 0, 40},
        "tritonAngle", "tritonEnergy");
    auto hKin_alfaAlfaTriton = dfAlfaTriton.Histo2D(
        {"hKin_alfaAlfaTriton", "Alfa Kinematics (alfa+triton);Theta [deg];Energy [MeV]", 100, 0, 180, 100, 0, 40},
        "alfaAngle", "alfaEnergy");

    auto hKin_tritonOnlyTriton = dfOnlyTriton.Histo2D(
        {"hKin_tritonOnlyTriton", "Triton Kinematics (only triton);Theta [deg];Energy [MeV]", 100, 0, 180, 100, 0, 40},
        "tritonAngle", "tritonEnergy");
    auto hKin_alfaOnlyAlfa = dfOnlyAlfa.Histo2D(
        {"hKin_alfaOnlyAlfa", "Alfa Kinematics (only alfa);Theta [deg];Energy [MeV]", 100, 0, 180, 100, 0, 40},
        "alfaAngle", "alfaEnergy");

    // Get Ex plot for every case
    auto hEx_All = df.Histo1D({"hEx_All", "Excitation Energy;E_{ex} [MeV];Counts", 100, 0, 16}, "Eex");
    auto hEx_AlfaTriton = dfAlfaTriton.Histo1D(
        {"hEx_AlfaTriton", "Excitation Energy (alfa+triton);E_{ex} [MeV];Counts", 100, 0, 16}, "Eex");
    auto hEx_OnlyTriton = dfOnlyTriton.Histo1D(
        {"hEx_OnlyTriton", "Excitation Energy (only triton);E_{ex} [MeV];Counts", 100, 0, 16}, "Eex");
    auto hEx_OnlyAlfa =
        dfOnlyAlfa.Histo1D({"hEx_OnlyAlfa", "Excitation Energy (only alfa);E_{ex} [MeV];Counts", 100, 0, 16}, "Eex");

    auto* cKinProton = new TCanvas("cKinProton", "Proton Kinematics", 1200, 800);
    cKinProton->Divide(2, 2);
    cKinProton->cd(1);
    hKin_protonAll->DrawClone("colz");
    cKinProton->cd(2);
    hKin_protonOnlyTriton->DrawClone("colz");
    cKinProton->cd(3);
    hKin_protonOnlyAlfa->DrawClone("colz");
    cKinProton->cd(4);
    hKin_protonAlfaTriton->DrawClone("colz");

    auto* cKinTritonAndAlfa = new TCanvas("cKinTritonAndAlfa", "Triton and Alfa Kinematics", 1200, 800);
    cKinTritonAndAlfa->Divide(3, 2);
    cKinTritonAndAlfa->cd(1);
    hKin_tritonAll->DrawClone("colz");
    cKinTritonAndAlfa->cd(2);
    hKin_tritonOnlyTriton->DrawClone("colz");
    cKinTritonAndAlfa->cd(3);
    hKin_tritonAlfaTriton->DrawClone("colz");
    cKinTritonAndAlfa->cd(4);
    hKin_alfaAll->DrawClone("colz");
    cKinTritonAndAlfa->cd(5);
    hKin_alfaOnlyAlfa->DrawClone("colz");
    cKinTritonAndAlfa->cd(6);
    hKin_alfaAlfaTriton->DrawClone("colz");

    auto* cEx = new TCanvas("cEx", "Excitation Energy", 1200, 800);
    cEx->Divide(2, 2);
    cEx->cd(1);
    hEx_All->DrawClone("hist");
    cEx->cd(2);
    hEx_OnlyTriton->DrawClone("hist");
    cEx->cd(3);
    hEx_OnlyAlfa->DrawClone("hist");
    cEx->cd(4);
    hEx_AlfaTriton->DrawClone("hist");
}