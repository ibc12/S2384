#include "ActDecayGenerator.h"
#include "ActKinematicGenerator.h"
#include "ActKinematics.h"
#include "ActParticle.h"
#include "ActSRIM.h"

#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TROOT.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// Alias corto para no repetir el tipo de nodo "tipo-borrado" de RDataFrame
using RNode = ROOT::RDF::RNode;

// Crea un TH2D (theta vs energia) para un df dado
ROOT::RDF::RResultPtr<TH2D> BookKin2D(RNode d, const std::string& name, const std::string& title,
                                      const std::string& thetaVar, const std::string& eVar, double eMax)
{
    return d.Histo2D({name.c_str(), title.c_str(), 100, 0, 180, 100, 0, eMax}, thetaVar, eVar);
}

// Crea un TH1D de Eex para un df dado
ROOT::RDF::RResultPtr<TH1D> BookEx1D(RNode d, const std::string& name, const std::string& title)
{
    return d.Histo1D({name.c_str(), title.c_str(), 100, 0, 16}, "Eex");
}

void analysisDecaySimulation()
{
    // Files for proton channel (light = 1H)
    // ROOT::RDataFrame df {"SimulationTTree",
    //                      "../Outputs/7Li/Decay/2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_decay_democratic4Body.root"};
    // ROOT::RDataFrame dfL1 {"LightMissTree",
    //                        "../Outputs/7Li/Decay/2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_decay_democratic4Body.root"};
    // Files for deuteron channel (light = 2H)
    ROOT::RDataFrame df {"SimulationTTree",
                         "../Outputs/7Li/Decay/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_decay_democratic3Body_dd.root"};
    ROOT::RDataFrame dfL1 {"LightMissTree",
                           "../Outputs/7Li/Decay/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_decay_democratic3Body_dd.root"};


    // ---- Selecciones sobre el arbol principal ----
    RNode dfTriton = df.Filter("tritonLayerCode != 0", "Select events where triton is detected");
    RNode dfAlfa = df.Filter("alfaLayerCode != 0", "Select events where alfa is detected");
    RNode dfAlfaTriton = df.Filter("gateAlfaTriton != 0", "Select events where both alfa and triton are detected");
    RNode dfOnlyTriton =
        df.Filter("tritonLayerCode != 0 && gateAlfaTriton == 0", "Select events where triton is detected but not alfa");
    RNode dfOnlyAlfa =
        df.Filter("alfaLayerCode != 0 && gateAlfaTriton == 0", "Select events where alfa is detected but not triton");

    // ---- Selecciones sobre LightMissTree (caso L1) ----
    RNode dfL1OnlyTriton = dfL1.Filter("tritonLayerCode != 0 && gateAlfaTriton == 0",
                                       "Select events where triton is detected but not alfa in L1");
    RNode dfL1OnlyAlfa = dfL1.Filter("alfaLayerCode != 0 && gateAlfaTriton == 0",
                                     "Select events where alfa is detected but not triton in L1");
    RNode dfL1AlfaTriton =
        dfL1.Filter("gateAlfaTriton != 0", "Select events where both alfa and triton are detected in L1");

    // Start by reading some SRIM tables to get the ranges of the particles in the gas
    auto* srim {new ActPhysics::SRIM};
    srim->ReadTable("alfaGas", "../../Calibrations/SRIM/4He_900mb_CF4_95-5.txt");
    srim->ReadTable("tritonGas", "../../Calibrations/SRIM/3H_900mb_CF4_95-5.txt");

    auto minRangeAlfa {srim->EvalRange("alfaGas", 6.0)};     // MeV
    auto minRangeTriton {srim->EvalRange("tritonGas", 4.0)}; // MeV

    std::cout << "Minimum range for alfa: " << minRangeAlfa << " mm" << std::endl;
    std::cout << "Minimum range for triton: " << minRangeTriton << " mm" << std::endl;

    // =========================================================
    // Caso principal (SimulationTTree)
    // =========================================================

    // Proton (light) kinematics
    auto hKin_lightAll =
        BookKin2D(df, "hKin_lightAll", "Light Kinematics;Theta [deg];Energy [MeV]", "theta3Lab", "EVertex", 20);
    auto hKin_lightOnlyTriton =
        BookKin2D(dfOnlyTriton, "hKin_lightOnlyTriton", "Light Kinematics (only triton);Theta [deg];Energy [MeV]",
                  "theta3Lab", "EVertex", 20);
    auto hKin_lightOnlyAlfa =
        BookKin2D(dfOnlyAlfa, "hKin_lightOnlyAlfa", "Light Kinematics (only alfa);Theta [deg];Energy [MeV]",
                  "theta3Lab", "EVertex", 20);
    auto hKin_lightAlfaTriton =
        BookKin2D(dfAlfaTriton, "hKin_lightAlfaTriton", "Light Kinematics (alfa+triton);Theta [deg];Energy [MeV]",
                  "theta3Lab", "EVertex", 20);

    // Triton / alfa kinematics
    auto hKin_tritonAll = BookKin2D(df, "hKin_tritonAll", "Triton Kinematics;Theta [deg];Energy [MeV]", "tritonAngle",
                                    "tritonEnergy", 40);
    auto hKin_alfaAll =
        BookKin2D(df, "hKin_alfaAll", "Alfa Kinematics;Theta [deg];Energy [MeV]", "alfaAngle", "alfaEnergy", 40);

    auto hKin_tritonAlfaTriton =
        BookKin2D(dfAlfaTriton, "hKin_tritonAlfaTriton", "Triton Kinematics (alfa+triton);Theta [deg];Energy [MeV]",
                  "tritonAngle", "tritonEnergy", 40);
    auto hKin_alfaAlfaTriton =
        BookKin2D(dfAlfaTriton, "hKin_alfaAlfaTriton", "Alfa Kinematics (alfa+triton);Theta [deg];Energy [MeV]",
                  "alfaAngle", "alfaEnergy", 40);

    auto hKin_tritonOnlyTriton =
        BookKin2D(dfOnlyTriton, "hKin_tritonOnlyTriton", "Triton Kinematics (only triton);Theta [deg];Energy [MeV]",
                  "tritonAngle", "tritonEnergy", 40);
    auto hKin_alfaOnlyAlfa =
        BookKin2D(dfOnlyAlfa, "hKin_alfaOnlyAlfa", "Alfa Kinematics (only alfa);Theta [deg];Energy [MeV]", "alfaAngle",
                  "alfaEnergy", 40);

    // Eex
    auto hEx_All = BookEx1D(df, "hEx_All", "Excitation Energy;E_{ex} [MeV];Counts");
    auto hEx_AlfaTriton =
        BookEx1D(dfAlfaTriton, "hEx_AlfaTriton", "Excitation Energy (alfa+triton);E_{ex} [MeV];Counts");
    auto hEx_OnlyTriton =
        BookEx1D(dfOnlyTriton, "hEx_OnlyTriton", "Excitation Energy (only triton);E_{ex} [MeV];Counts");
    auto hEx_OnlyAlfa = BookEx1D(dfOnlyAlfa, "hEx_OnlyAlfa", "Excitation Energy (only alfa);E_{ex} [MeV];Counts");

    // =========================================================
    // Caso L1 (LightMissTree)
    // =========================================================

    // "Light" (proton perdido) kinematics -> aqui la energia es T3Lab, no EVertex
    auto hKin_lightAll_L1 =
        BookKin2D(dfL1, "hKin_lightAll_L1", "Light Kinematics L1;Theta [deg];Energy [MeV]", "theta3Lab", "T3Lab", 20);
    auto hKin_lightOnlyTriton_L1 =
        BookKin2D(dfL1OnlyTriton, "hKin_lightOnlyTriton_L1",
                  "Light Kinematics L1 (only triton);Theta [deg];Energy [MeV]", "theta3Lab", "T3Lab", 20);
    auto hKin_lightOnlyAlfa_L1 =
        BookKin2D(dfL1OnlyAlfa, "hKin_lightOnlyAlfa_L1", "Light Kinematics L1 (only alfa);Theta [deg];Energy [MeV]",
                  "theta3Lab", "T3Lab", 20);
    auto hKin_lightAlfaTriton_L1 =
        BookKin2D(dfL1AlfaTriton, "hKin_lightAlfaTriton_L1",
                  "Light Kinematics L1 (alfa+triton);Theta [deg];Energy [MeV]", "theta3Lab", "T3Lab", 20);

    // Triton / alfa kinematics en L1
    auto hKin_tritonAll_L1 = BookKin2D(dfL1, "hKin_tritonAll_L1", "Triton Kinematics L1;Theta [deg];Energy [MeV]",
                                       "tritonAngle", "tritonEnergy", 40);
    auto hKin_alfaAll_L1 = BookKin2D(dfL1, "hKin_alfaAll_L1", "Alfa Kinematics L1;Theta [deg];Energy [MeV]",
                                     "alfaAngle", "alfaEnergy", 40);

    auto hKin_tritonAlfaTriton_L1 =
        BookKin2D(dfL1AlfaTriton, "hKin_tritonAlfaTriton_L1",
                  "Triton Kinematics L1 (alfa+triton);Theta [deg];Energy [MeV]", "tritonAngle", "tritonEnergy", 40);
    auto hKin_alfaAlfaTriton_L1 =
        BookKin2D(dfL1AlfaTriton, "hKin_alfaAlfaTriton_L1", "Alfa Kinematics L1 (alfa+triton);Theta [deg];Energy [MeV]",
                  "alfaAngle", "alfaEnergy", 40);

    auto hKin_tritonOnlyTriton_L1 =
        BookKin2D(dfL1OnlyTriton, "hKin_tritonOnlyTriton_L1",
                  "Triton Kinematics L1 (only triton);Theta [deg];Energy [MeV]", "tritonAngle", "tritonEnergy", 40);
    auto hKin_alfaOnlyAlfa_L1 =
        BookKin2D(dfL1OnlyAlfa, "hKin_alfaOnlyAlfa_L1", "Alfa Kinematics L1 (only alfa);Theta [deg];Energy [MeV]",
                  "alfaAngle", "alfaEnergy", 40);

    // Eex en L1
    auto hEx_All_L1 = BookEx1D(dfL1, "hEx_All_L1", "Excitation Energy L1;E_{ex} [MeV];Counts");
    auto hEx_AlfaTriton_L1 =
        BookEx1D(dfL1AlfaTriton, "hEx_AlfaTriton_L1", "Excitation Energy L1 (alfa+triton);E_{ex} [MeV];Counts");
    auto hEx_OnlyTriton_L1 =
        BookEx1D(dfL1OnlyTriton, "hEx_OnlyTriton_L1", "Excitation Energy L1 (only triton);E_{ex} [MeV];Counts");
    auto hEx_OnlyAlfa_L1 =
        BookEx1D(dfL1OnlyAlfa, "hEx_OnlyAlfa_L1", "Excitation Energy L1 (only alfa);E_{ex} [MeV];Counts");

    // =========================================================
    // Canvases: caso principal
    // =========================================================
    auto* cKinLight = new TCanvas("cKinLight", "Light Kinematics", 1200, 800);
    cKinLight->Divide(2, 2);
    cKinLight->cd(1);
    hKin_lightAll->DrawClone("colz");
    cKinLight->cd(2);
    hKin_lightOnlyTriton->DrawClone("colz");
    cKinLight->cd(3);
    hKin_lightOnlyAlfa->DrawClone("colz");
    cKinLight->cd(4);
    hKin_lightAlfaTriton->DrawClone("colz");

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

    // =========================================================
    // Canvases: caso L1
    // =========================================================
    auto* cKinLight_L1 = new TCanvas("cKinLight_L1", "Light Kinematics (L1)", 1200, 800);
    cKinLight_L1->Divide(2, 2);
    cKinLight_L1->cd(1);
    hKin_lightAll_L1->DrawClone("colz");
    cKinLight_L1->cd(2);
    hKin_lightOnlyTriton_L1->DrawClone("colz");
    cKinLight_L1->cd(3);
    hKin_lightOnlyAlfa_L1->DrawClone("colz");
    cKinLight_L1->cd(4);
    hKin_lightAlfaTriton_L1->DrawClone("colz");

    auto* cKinTritonAndAlfa_L1 = new TCanvas("cKinTritonAndAlfa_L1", "Triton and Alfa Kinematics (L1)", 1200, 800);
    cKinTritonAndAlfa_L1->Divide(3, 2);
    cKinTritonAndAlfa_L1->cd(1);
    hKin_tritonAll_L1->DrawClone("colz");
    cKinTritonAndAlfa_L1->cd(2);
    hKin_tritonOnlyTriton_L1->DrawClone("colz");
    cKinTritonAndAlfa_L1->cd(3);
    hKin_tritonAlfaTriton_L1->DrawClone("colz");
    cKinTritonAndAlfa_L1->cd(4);
    hKin_alfaAll_L1->DrawClone("colz");
    cKinTritonAndAlfa_L1->cd(5);
    hKin_alfaOnlyAlfa_L1->DrawClone("colz");
    cKinTritonAndAlfa_L1->cd(6);
    hKin_alfaAlfaTriton_L1->DrawClone("colz");

    auto* cEx_L1 = new TCanvas("cEx_L1", "Excitation Energy (L1)", 1200, 800);
    cEx_L1->Divide(2, 2);
    cEx_L1->cd(1);
    hEx_All_L1->DrawClone("hist");
    cEx_L1->cd(2);
    hEx_OnlyTriton_L1->DrawClone("hist");
    cEx_L1->cd(3);
    hEx_OnlyAlfa_L1->DrawClone("hist");
    cEx_L1->cd(4);
    hEx_AlfaTriton_L1->DrawClone("hist");
}