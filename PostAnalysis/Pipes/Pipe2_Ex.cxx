#ifndef Pipe2_Ex_cxx
#define Pipe2_Ex_cxx

#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActParticle.h"
#include "ActSRIM.h"

#include "ROOT/RDataFrame.hxx"
#include "Rtypes.h"

#include "TCanvas.h"
#include "TMath.h"
#include "TROOT.h"
#include "TString.h"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "../HistConfig.h"

void Pipe2_Ex(const std::string& beam, const std::string& target, const std::string& light)
{
    // Read data
    auto filename {TString::Format("./Outputs/tree_pid_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str())};
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"PID_Tree", filename};

    // Init SRIM
    auto* srim {new ActPhysics::SRIM};
    // Correct SRIM names
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
    srim->ReadTable(light, TString::Format("../Calibrations/SRIM/%s_900mb_CF4_95-5.txt", srimName.c_str()).Data());
    srim->ReadTable(beam, TString::Format("../Calibrations/SRIM/%s_900mb_CF4_95-5.txt", beam.c_str()).Data());
    // Build energy at vertex
    auto dfVertex = df.Define("EVertex",
                              [&](const ActRoot::MergerData& d)
                              {
                                  double ret {};
                                  if(d.fLight.IsFilled())
                                      ret = srim->EvalInitialEnergy(light, d.fLight.fEs.front(), d.fLight.fTL);
                                  else // L1 trigger
                                      ret = srim->EvalEnergy(light, d.fLight.fTL);
                                  return ret;
                              },
                              {"MergerData"});

    // Init particles
    ActPhysics::Particle pb {beam};
    ActPhysics::Particle pt {target};
    ActPhysics::Particle pl {light};

    // // Filter on heavy particle hit in the telescope
    auto def {dfVertex};
    // auto def = dfVertex.Filter([](const ActRoot::MergerData &m)
    //                            { if(!m.fHeavy.fLayers.empty() && m.fHeavy.fLayers.front() == "f2")
    //                            {
    //                                 if(m.fHeavy.fEs[0] > 9.5)
    //                                 {
    //                                     return true;
    //                                 }
    //                                 else
    //                                     return false;
    //                            }
    //
    //                            else
    //                                 return false; }, {"MergerData"});
    //
    // Build beam energy
    def = def.Define("EBeam", [&](const ActRoot::MergerData& d)
                     { return srim->Slow(beam, 7.5 * pb.GetAMU(), d.fRP.X()); }, {"MergerData"});

    ActPhysics::Kinematics kin {pb, pt, pl, 7.5 * pb.GetAMU()};
    // Vector of kinematics as one object is needed per
    // processing slot (since we are changing EBeam in each entry)
    std::vector<ActPhysics::Kinematics> vkins {def.GetNSlots()};
    for(auto& k : vkins)
        k = kin;
    def =
        def.DefineSlot("Ex",
                       [&](unsigned int slot, const ActRoot::MergerData& d, double EVertex, double EBeam)
                       {
                           vkins[slot].SetBeamEnergy(EBeam);
                           return vkins[slot].ReconstructExcitationEnergy(EVertex, (d.fThetaLight) * TMath::DegToRad());
                       },
                       {"MergerData", "EVertex", "EBeam"});
    def =
        def.DefineSlot("ThetaCM",
                       [&](unsigned int slot, const ActRoot::MergerData& d, double EVertex, double EBeam)
                       {
                           vkins[slot].SetBeamEnergy(EBeam);
                           return vkins[slot].ReconstructTheta3CMFromLab(EVertex, (d.fThetaLight) * TMath::DegToRad()) *
                                  TMath::RadToDeg();
                       },
                       {"MergerData", "EVertex", "EBeam"});


    // Book new histograms
    auto hKin {def.Histo2D(HistConfig::KinEl, "fThetaLight", "EVertex")};

    auto hKinCM {def.Histo2D(HistConfig::KinCM, "ThetaCM", "EVertex")};

    auto hEBeam {def.Histo1D("EBeam")};

    auto hExSil {def.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == true; }, {"MergerData"})
                     .Histo1D(HistConfig::Ex, "Ex")};
    hExSil->SetTitle("Ex with silicons");
    auto hExL1 {def.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"})
                    .Histo1D(HistConfig::Ex, "Ex")};
    hExL1->SetTitle("Ex with L1");

    auto hTheta {def.Histo1D("fThetaLight")};

    auto hThetaBeam {def.Histo2D(HistConfig::ThetaBeam, "fRP.fCoordinates.fX", "fThetaBeam")};

    auto hRP {def.Histo2D(HistConfig::RP, "fRP.fCoordinates.fX", "fRP.fCoordinates.fY")};

    auto hThetaCMLab {def.Histo2D(HistConfig::ThetaCMLab, "fThetaLight", "ThetaCM")};

    // Ex dependences
    auto hExThetaCM {def.Histo2D(HistConfig::ExThetaCM, "ThetaCM", "Ex")};
    auto hExThetaLab {def.Histo2D(HistConfig::ExThetaLab, "fThetaLight", "Ex")};
    auto hExRP {def.Histo2D(HistConfig::ExRPx, "fRP.fCoordinates.fX", "Ex")};

    // Heavy histograms
    auto hThetaHLLab {def.Histo2D(HistConfig::ChangeTitle(HistConfig::ThetaHeavyLight, "Lab correlations"),
                                  "fThetaLight", "fThetaHeavy")};
    // Save!
    auto outfile {TString::Format("./Outputs/tree_ex_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str())};
    def.Snapshot("Final_Tree", outfile);
    std::cout << "Saving Final_Tree in " << outfile << '\n';
    // Save Ex histos on file
    auto file {std::make_shared<TFile>(outfile.Data(), "update")};
    hExSil->Write("hExSil");
    hExL1->Write("hExL1");
    file->Close();


    auto* c22 {new TCanvas("c22", "Pipe2 canvas 2")};
    c22->DivideSquare(4);
    c22->cd(1);
    hTheta->DrawClone();
    c22->cd(2);
    hThetaBeam->DrawClone("colz");
    c22->cd(3);
    hRP->DrawClone("colz");
    c22->cd(4);
    hEBeam->DrawClone();

    auto* c21 {new TCanvas("c21", "Pipe2 canvas 1")};
    c21->DivideSquare(6);
    c21->cd(1);
    hKin->DrawClone("colz");
    auto* theo {kin.GetKinematicLine3()};
    theo->Draw("same");
    c21->cd(2);
    hExSil->DrawClone();
    c21->cd(3);
    hKinCM->DrawClone("colz");
    c21->cd(4);
    hExL1->DrawClone();
    // hExThetaLab->DrawClone("colz");
    c21->cd(5);
    hExThetaCM->DrawClone("colz");
    c21->cd(6);
    hExRP->DrawClone("colz");

    auto* c23 {new TCanvas {"c23", "Pipe2 canvas 3"}};
    c23->DivideSquare(4);
    c23->cd(1);
    hThetaHLLab->DrawClone("colz");
    c23->cd(2);
    hThetaCMLab->DrawClone("colz");
    c23->cd(3);
}
#endif
