#ifndef triumf_all_cxx
#define triumf_all_cxx
#include "ActColors.h"
#include "ActCrossSection.h"
#include "ActDecayGenerator.h"
#include "ActKinematicGenerator.h"
#include "ActKinematics.h"
#include "ActLine.h"
#include "ActParticle.h"
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
#include "TMath.h"
#include "TRandom.h"
#include "TString.h"
#include "TTree.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>

#include "../PostAnalysis/HistConfig.h"

using XYZPoint = ROOT::Math::XYZPoint;
using XYZVector = ROOT::Math::XYZVector;

XYZPoint SampleVertex(ActRoot::TPCParameters* tpc)
{
    // Define sigmas along Y and Z
    double sigmaY {4};
    double sigmaZ {4};
    auto y {gRandom->Gaus(tpc->Y() / 2, sigmaY)};
    auto z {gRandom->Gaus(135., sigmaZ)}; // Actar has the beam entrance 135 mm from the bottom of the field cage.
    auto x {gRandom->Uniform() * tpc->X()};
    return {x, y, z};
}

std::pair<double, double> SampleCM()
{
    auto theta {TMath::ACos(gRandom->Uniform(-1, 1))};
    auto phi {gRandom->Uniform(0, TMath::TwoPi())};
    return {theta, phi};
}

void ApplyNaN(double& e, double t = 0, const std::string& comment = "stopped")
{
    if(e <= t)
        e = std::nan(comment.c_str());
}

void ApplyThetaRes(double& theta)
{
    double sigma {0.95 / 2.355}; // FWHM to sigma
    theta = gRandom->Gaus(theta, sigma * TMath::DegToRad());
}

double RandomizeBeamEnergy(double Tini, double sigma)
{
    return gRandom->Gaus(Tini, sigma);
}

std::map<std::string, double> LoadEfficiencies(const std::string& filename)
{
    std::map<std::string, double> efficiencies;
    std::ifstream fin(filename);
    std::string key;
    double value;
    while(fin >> key >> value)
    {
        efficiencies[key] = value;
    }
    fin.close();
    return efficiencies;
}

// Some silicons malfunctioned in some of the experimental runs, this takes into account the ammount of time they were
// not working
bool AcceptHit(std::map<std::string, double> efficiencies, const std::string& layer, int detID)
{
    TString key = Form("%s_%d", layer.c_str(), detID);
    auto it = efficiencies.find(key.Data());
    if(it == efficiencies.end())
        return true; // si no está definido → aceptar
    double eff = it->second;
    return (gRandom->Rndm() < eff);
}

// Get XS files depending on the reaction in place
bool GetXS(const std::string& target, const std::string& light, const std::string& beam, double Ex,
           ActSim::CrossSection* xs)
{
    bool isThereXS {};
    if(target == "2H" && light == "1H" && beam == "11Li")
    {
        isThereXS = true;
        if(Ex == 0.)
        {
            TString data_to_read {TString::Format("./Inputs/xs/%s/dp/angs12nospin.dat", beam.c_str())};
            xs->ReadFile(data_to_read.Data());
            std::cout << "Total xs: " << xs->GetTotalXSmbarn() << std::endl;
        }
        else if(Ex == 0.130)
        {
            TString data_to_read {TString::Format("./Inputs/xs/%s/dp/angp12nospin.dat", beam.c_str())};
            xs->ReadFile(data_to_read.Data());
            std::cout << "Total xs: " << xs->GetTotalXSmbarn() << std::endl;
        }
        else if(Ex == 0.435)
        {
            TString data_to_read {TString::Format("./Inputs/xs/%s/dp/angp32nospin.dat", beam.c_str())};
            xs->ReadFile(data_to_read.Data());
            std::cout << "Total xs: " << xs->GetTotalXSmbarn() << std::endl;
        }
        else if(Ex == 2.)
        {
            TString data_to_read {TString::Format("./Inputs/xs/%s/dp/angd52nospin.dat", beam.c_str())};
            xs->ReadFile(data_to_read.Data());
            std::cout << "Total xs: " << xs->GetTotalXSmbarn() << std::endl;
        }
        else if(Ex == 5.)
        {
            TString data_to_read {TString::Format("./Inputs/xs/%s/dp/angd52nospin.dat", beam.c_str())};
            xs->ReadFile(data_to_read.Data());
            std::cout << "Total xs: " << xs->GetTotalXSmbarn() << std::endl;
        }
    }
    else if(target == "2H" && light == "2H" && beam == "11Li")
    {
        if(Ex == 0.)
        {
            isThereXS = true;
            TString data_to_read {TString::Format("./Inputs/xs/%s/dd/elastic.dat", beam.c_str())};
            xs->ReadFile(data_to_read.Data());
            std::cout << "Total xs: " << xs->GetTotalXSmbarn() << std::endl;
        }
    }
    else if(target == "2H" && light == "2H" && beam == "7Li")
    {
        if(Ex == 0.)
        {
            isThereXS = true;
            TString data_to_read {TString::Format("./Inputs/xs/%s/dd/elastic.dat", beam.c_str())};
            xs->ReadFile(data_to_read.Data());
            std::cout << "Total xs: " << xs->GetTotalXSmbarn() << std::endl;
        }
    }
    return isThereXS;
}

void CheckL1Acceptance(XYZVector direction, XYZPoint vertex, XYZPoint finalPointgas, double minPads,
                       double halfWidthExclusionZone)
{
    int a = 1;
}

void do_all_simus(const std::string& beam, const std::string& target, const std::string& light,
                  const std::string& heavy, int neutronPS, int protonPS, double Tbeam, double Ex, bool inspect)
{
    // Set number of iterations
    auto niter {static_cast<int>(1e7)};
    gRandom->SetSeed(0);
    // Initialize detectors
    // TPC
    ActRoot::TPCParameters tpc {"Actar"};
    std::cout << "TPC: " << tpc.X() << " " << tpc.Y() << " " << tpc.Z() << '\n';
    // Silicons
    auto* sils {new ActPhysics::SilSpecs};
    std::string silConfig("silspecs"); // no front silicons, only lateral ones
    sils->ReadFile("../configs/" + silConfig + ".conf");
    sils->Print();
    const double sigmaSil {0.100 / 2.355}; // Si resolution for laterals, around 100 keV FWHM
    auto silRes = std::make_unique<TF1>(
        "silRes", [=](double* x, double* p) { return sigmaSil * TMath::Sqrt(x[0] / 5.5); }, 0.0, 100.0, 1);
    std::vector<std::string> silLayers {"f0", "l0", "r0"};
    std::vector<std::string> AllsilLayers {"f0", "f1", "f2", "f3", "l0", "r0"};
    // We have to centre the silicons with the beam input
    // In real life beam window is not at Z / 2
    for(auto& [name, layer] : sils->GetLayers())
    {
        if(name == "f0" || name == "f1")
            layer.MoveZTo(75, {3});
        if(name == "f2" || name == "f3")
            layer.MoveZTo(125, {0});
        if(name == "l0" || name == "r0") // beam went at the height of the half of the second highest silicon
            layer.MoveZTo(135, {7});
    }
    // Silicon malfunction txt
    std::string silEfficienciesPath {"./Inputs/Efficiencies/silicon_efficiencies.txt"};
    std::map<std::string, double> silEfficiencies {LoadEfficiencies(silEfficienciesPath)};

    std::cout << "Sils Z centred at : " << tpc.Z() / 2 << " mm" << '\n';
    sils->DrawGeo();
    // This means: make the Z of silicons {5,6,...} be that zOfBeam.
    // shift the others accordingly
    // sils->DrawGeo();

    // Sigmas
    const double sigmaPercentBeam {0.008};
    // Flags for resolution
    bool RestOfBeamLine {true}; // If true enables CFA and mylar of entrance
    bool exResolution {true};

    // SRIM
    auto* srim {new ActPhysics::SRIM};
    std::string path {"../Calibrations/SRIM/"};
    std::string gas {"900mb_CF4_95-5"};
    std::string CFAgas {"6mb_butane"};
    std::string Mylar {"Mylar"};
    std::string silicon {"silicon"};
    srim->ReadTable("beam", path + beam + "_" + gas + ".txt");
    // srim->ReadTable("beamCFA", path + beam + "_" + CFAgas + ".txt");
    srim->ReadTable("beamMylar", path + beam + "_" + Mylar + ".txt");
    srim->ReadTable("light", path + light + "_" + gas + ".txt");
    srim->ReadTable("heavy", path + heavy + "_" + gas + ".txt");
    srim->ReadTable("lightInSil", path + light + "_" + silicon + ".txt");
    srim->ReadTable("heavyInSil", path + heavy + "_" + silicon + ".txt");
    // srim->SetStragglingLISE("heavyInSil", "../Calibrations/LISE files/" + heavy + "_silicon" + ".txt");
    // srim->SetStragglingLISE("heavy", "../Calibrations/LISE files/" + heavy + "_gas_95-5" + ".txt");
    // srim->SetStragglingLISE("lightInSil", "../Calibrations/LISE files/" + light + "_silicon" + ".txt");
    // srim->SetStragglingLISE("light", "../Calibrations/LISE files/" + light + "_gas_95-5" + ".txt");
    // srim->SetStragglingLISE("beamMylar", "../Calibrations/LISE files/" + beam + "_Mylar" + ".txt");
    // srim->SetStragglingLISE("beamCFA", "../Calibrations/LISE files/" + beam + "_gasCFA" + ".txt");
    // srim->SetStragglingLISE("beam", "../Calibrations/LISE files/" + beam + "_gas_95-5" + ".txt");

    // Kinematics
    auto* kinTheo {new ActPhysics::Kinematics {beam, target, light, heavy, Tbeam, Ex}};
    auto* kin {new ActPhysics::Kinematics {beam, target, light, heavy, Tbeam, Ex}};
    auto* kinGen {new ActSim::KinematicGenerator {beam, target, light, heavy, protonPS, neutronPS}};

    // cross section sampler
    auto* xs {new ActSim::CrossSection()};
    bool isThereXS {};
    if(neutronPS == 0 && protonPS == 0)
    {
        isThereXS = GetXS(target, light, beam, Ex, xs);
    }
    double alpha {1.};
    double NLi11 {9.87839e8};              // Counts in CFA trigger corrected by CFA/F2 ratio
    double NLi7 {2.08039e8};               // Counts in CFA trigger corrected by CFA/F2 ratio
    double Nd {4.6688e19 * 25.6 * 0.8877}; // atom density, 25,6 cm long, 88.77% d2
    if(isThereXS)
    {
        if(beam == "7Li")
        {
            alpha = (NLi7 * Nd * xs->GetTotalXScm2()) / niter;
        }
        else if(beam == "11Li")
        {
            alpha = (NLi11 * Nd * xs->GetTotalXScm2()) / niter;
        }
        std::cout << "Alpha: " << alpha << std::endl;
    }
    // Declare histograms
    // kinematics and angles
    auto hKin {HistConfig::Kin.GetHistogram()};
    auto hTheta3CM {HistConfig::ThetaCM.GetHistogram()};
    auto hThetaCMAll {HistConfig::ThetaCM.GetHistogram()};
    hThetaCMAll->SetTitle("Theta3CM all;#theta_{CM} [#circ];Counts");
    auto hThetaLabAll {HistConfig::ThetaCM.GetHistogram()};
    hThetaLabAll->SetTitle("Theta3Lab all;#theta_{Lab} [#circ];Counts");
    auto hTheta3Lab {HistConfig::ThetaCM.GetHistogram()};
    hTheta3Lab->SetTitle("Theta3Lab;#theta_{Lab} [#circ];Counts");
    auto hPhiAll {HistConfig::PhiCM.GetHistogram()};
    hPhiAll->SetTitle("Phi3CM all;#phi_{CM} [#circ];Counts");
    auto hPhi3CM {HistConfig::PhiCM.GetHistogram()};
    hPhi3CM->SetTitle("Phi3CM;#phi_{CM} [#circ];Counts");
    // Silicon hits
    auto hSPf0 {HistConfig::SP.GetHistogram()};
    hSPf0->SetTitle("SP for f0");
    auto hSPl0 {HistConfig::SP.GetHistogram()};
    hSPl0->SetTitle("SP for l0");
    auto hSPr0 {HistConfig::SP.GetHistogram()};
    hSPr0->SetTitle("SP for r0");
    // Reconstructed histos
    auto hEx {HistConfig::Ex.GetHistogram()};
    auto hKinRec {HistConfig::Kin.GetHistogram()};
    hKinRec->SetTitle("Reconstructed Kinetic Energy;T_{3} [MeV];Counts");
    auto hRP_X {HistConfig::RPx.GetHistogram()};
    hRP_X->SetTitle("Reconstructed RP;X [mm];Counts");

    // File to save data
    TString fileName {TString::Format("./Outputs/%s/%s_%s_TRIUMF_Eex_%.3f_nPS_%d_pPS_%d.root", beam.c_str(),
                                      target.c_str(), light.c_str(), Ex, neutronPS, protonPS)};
    auto* outFile {new TFile(fileName, "recreate")};
    auto* outTree {new TTree("SimulationTTree", "A TTree containing only our Eex obtained by simulation")};
    double theta3CM_tree {};
    outTree->Branch("theta3CM", &theta3CM_tree);
    double Eex_tree {};
    outTree->Branch("Eex", &Eex_tree);
    double EVertex_tree {};
    outTree->Branch("EVertex", &EVertex_tree);
    double theta3Lab_tree {};
    outTree->Branch("theta3Lab", &theta3Lab_tree);
    double phi3CM_tree {};
    outTree->Branch("phi3CM", &phi3CM_tree);
    // Set Random Ex if needed (no xs available, so will be uniform distributed)
    if(neutronPS == 2)
    {
        Ex = (1.26642 + 0.36928) / 2;
    }
    // RUN!
    // print fancy info (dont pay attention to this)
    std::cout << BOLDMAGENTA << "Running for Ex = " << Ex << " MeV" << RESET << '\n';
    std::cout << BOLDGREEN;
    const int percentPrint {5};
    int step {niter / (100 / percentPrint)};
    int nextPrint {step};
    int percent {};
    for(int it = 0; it < niter; it++)
    {
        // Print progress
        if(it >= nextPrint)
        {
            percent = 100 * it / niter;
            std::cout << "\r" << std::string(percent / percentPrint, '|') << percent << "%";
            std::cout.flush();
            nextPrint += step;
        }
        // Sample vertex position
        auto vertex {SampleVertex(&tpc)};

        // Randomize (if needed) Ex in a BW distribution
        double randEx = Ex;
        if(exResolution && isThereXS && light == "1H" && beam == "11Li")
        {
            if(Ex == 0)
            {
                randEx = gRandom->BreitWigner(Ex, 0.1);
            }
            else if(Ex == 0.130)
            {
                randEx = gRandom->BreitWigner(Ex, 0.015);
            }
            else if(Ex == 0.435)
            {
                randEx = gRandom->BreitWigner(Ex, 0.08);
            }
            else if(Ex == 2)
            {
                randEx = gRandom->BreitWigner(Ex, 0.08);
            }
            else if(Ex == 5)
            {
                randEx = gRandom->BreitWigner(Ex, 0.08);
            }
        }
        // Randomize beam energy, slow beam with straggling and check if reaction can happen
        auto TbeamRand = RandomizeBeamEnergy(Tbeam, sigmaPercentBeam * Tbeam);
        if(RestOfBeamLine)
        {
            // TbeamRand = srim->SlowWithStraggling("beamCFA", TbeamRand, 19);       // Gas CFA
            TbeamRand = srim->SlowWithStraggling("beamMylar", TbeamRand, 0.0168);  // All mylar
            TbeamRand = srim->SlowWithStraggling("beam", TbeamRand, 60);          // Gas before pad plane
        }
        auto TbeamCorr {srim->SlowWithStraggling("beam", TbeamRand, vertex.X())};
        // Initialize variables for both methods, kinGen and kin
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
        // Sample kinematics, diferent method depending on existance of xs and particles in ps
        if(isThereXS)
        {
            auto beamThreshold {ActPhysics::Kinematics(beam, target, light, heavy, -1, randEx).GetT1Thresh()};
            if(std::isnan(TbeamCorr) || TbeamCorr < beamThreshold)
            {
                continue;
            }
            kin->SetBeamEnergyAndEx(TbeamCorr, randEx);
            // Sample angle with xs
            while(theta3CMBefore < 0)
            {
                // theta3CMBefore = xs->SampleCDF(gRandom->Uniform());
                theta3CMBefore = xs->SampleHist();
                // std::cout << theta3CMBefore << std::endl;
            } // sample in deg
            phi3CM = gRandom->Uniform(0, 2 * TMath::Pi());
            kin->ComputeRecoilKinematics(theta3CMBefore * TMath::DegToRad(), phi3CM);
            // Get Lab kinematics
            T3Lab = kin->GetT3Lab();
            phi3Lab = kin->GetPhi3Lab();
            theta3Lab = kin->GetTheta3Lab();
            // Save without resolution
            theta3LabSampled = theta3Lab;
            // Apply angle resolution
            ApplyThetaRes(theta3Lab);
            theta3CM = kin->ReconstructTheta3CMFromLab(TbeamCorr, theta3Lab);

            // Heavy
            theta4Lab = kin->GetTheta4Lab();
            phi4Lab = kin->GetPhi4Lab();
            T4Lab = kin->GetT4Lab();
        }
        else
        {
            // Sample kinematics generator
            kinGen->SetBeamAndExEnergies(TbeamCorr, randEx);
            weight = kinGen->Generate();
            if(neutronPS == 0 && protonPS == 0)
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
            // Save without resolution
            theta3LabSampled = theta3Lab;
            // Apply angle resolution
            ApplyThetaRes(theta3Lab);

            theta3CMBefore = kin->ReconstructTheta3CMFromLab(T3Lab, theta3LabSampled) *
                             TMath::RadToDeg(); // this is in deg, because of xs sampling in other case
            theta3CM = kin->ReconstructTheta3CMFromLab(T3Lab, theta3Lab);
            phi3CM = phi3Lab;

            // Heavy
            theta4Lab = LorenztVector4->Theta();
            phi4Lab = LorenztVector4->Phi();
            T4Lab = LorenztVector4->E() - LorenztVector4->M();
        }
        // Fill kinematics and angles
        hKin->Fill(theta3LabSampled * TMath::RadToDeg(), T3Lab);
        hThetaCMAll->Fill(theta3CMBefore);
        hThetaLabAll->Fill(theta3Lab * TMath::RadToDeg());
        hPhiAll->Fill(phi3Lab * TMath::RadToDeg());
        // Extract direction
        XYZVector direction {TMath::Cos(theta3Lab), TMath::Sin(theta3Lab) * TMath::Sin(phi3Lab),
                             TMath::Sin(theta3Lab) * TMath::Cos(phi3Lab)};
        // Threshold L1, particles that stop in actar. Check before doing the continues
        double rangeInGas {srim->EvalRange("light", T3Lab)};
        ROOT::Math::XYZPoint finalPointGas {vertex + rangeInGas * direction.Unit()};
        if(0 <= finalPointGas.X() && finalPointGas.X() <= 256 && 0 <= finalPointGas.Y() && finalPointGas.Y() <= 256 &&
           0 <= finalPointGas.Z() && finalPointGas.Z() <= 256)
        {
        }
        // How to check whether tracks would read the silicons with new class:
        int silIndex0 = -1;
        ROOT::Math::XYZPoint silPoint0;
        std::string layer0;
        for(auto layer : silLayers)
        {
            std::tie(silIndex0, silPoint0) = sils->FindSPInLayer(layer, vertex, direction);
            if(silIndex0 != -1)
            {
                layer0 = layer;
                break;
            }
        }
        if(layer0 == "f0")
            continue; // for this simulation we don't consider front silicons
        // Check if corresponds to a hit when the detector was off or on
        if(!AcceptHit(silEfficiencies, layer0, silIndex0))
        {
            continue; // if not accepted, go to next iteration
        }
        if(silIndex0 == -1)
        {
            continue;
        }
        // Slow down light in gas
        auto T3AtSil {srim->SlowWithStraggling("light", T3Lab, (silPoint0 - vertex).R())};
        // Check if stopped
        ApplyNaN(T3AtSil);
        if(std::isnan(T3AtSil))
        {
            continue;
        }
        // Slow down in silicon
        auto normal {sils->GetLayer(layer0).GetNormal()};
        auto angleWithNormal {TMath::ACos(direction.Unit().Dot(normal.Unit()))};
        auto T3AfterSil0 {srim->SlowWithStraggling("lightInSil", T3AtSil,
                                                   sils->GetLayer(layer0).GetUnit().GetThickness(), angleWithNormal)};
        auto eLoss0preSilRes {T3AtSil - T3AfterSil0};
        auto eLoss0 {gRandom->Gaus(eLoss0preSilRes, silRes->Eval(eLoss0preSilRes))}; // after silicon resolution
        ApplyNaN(eLoss0, sils->GetLayer(layer0).GetThresholds().at(silIndex0));
        int count = 0;
        if(std::isnan(eLoss0))
        {
            continue;
        }

        // Apply 2nd layer of silicons
        double T3AfterInterGas {};
        int silIndex1 {};
        ROOT::Math::XYZPoint silPoint1 {};
        double eLoss1 {};
        double T3AfterSil1 {-1};
        if(T3AfterSil0 > 0. && layer0 == "f0")
        {
            std::tie(silIndex1, silPoint1) = sils->FindSPInLayer("f1", vertex, direction);
            if(silIndex1 == -1)
            {
            } // If a silicon is not reached, don't continue with punchthough calculation
            else
            {
                T3AfterInterGas = {srim->SlowWithStraggling("light", T3AfterSil0, (silPoint0 - silPoint1).R())};
                if(T3AfterInterGas == 0)
                {
                } // If slow in gas don't continue with calculation
                else
                {
                    T3AfterSil1 = srim->SlowWithStraggling(
                        "lightInSil", T3AfterInterGas, sils->GetLayer("f1").GetUnit().GetThickness(), angleWithNormal);
                    auto eLoss1preSilRes {T3AfterInterGas - T3AfterSil1};
                    eLoss1 = gRandom->Gaus(eLoss1preSilRes, silRes->Eval(eLoss1preSilRes)); // after silicon resolution
                    ApplyNaN(eLoss1, sils->GetLayer("f1").GetThresholds().at(silIndex1));
                    if(std::isnan(eLoss1))
                        eLoss1 = 0;
                }
            }
        }
        // Reconstruct!
        bool isOk {(T3AfterSil0 == 0 || T3AfterSil1 == 0)}; // no punchthrouhg
        if(isOk)
        {
            // Assuming no punchthrough!
            double T3Rec {};
            if(eLoss1 == 0)
            {
                T3Rec = srim->EvalInitialEnergy("light", eLoss0, (silPoint0 - vertex).R());
            }
            else
            {
                auto T3Rec0 {srim->EvalInitialEnergy("light", eLoss0, (silPoint0 - vertex).R())};
                // Reconstruction ob T3 with 2 silicon layers
                auto T3Rec1 {srim->EvalInitialEnergy("light", eLoss1, (silPoint1 - silPoint0).R())};
                T3Rec = srim->EvalInitialEnergy("light", eLoss0 + T3Rec1, (silPoint0 - vertex).R());
            }
            auto ExRec {kin->ReconstructExcitationEnergy(T3Rec, theta3Lab)};
            // Fill
            hKinRec->Fill(theta3Lab * TMath::RadToDeg(), T3Rec); // after reconstruction
            hEx->Fill(ExRec, weight); // To get real counts weigth * alpha
            hRP_X->Fill(vertex.X());
            if(layer0 == "f0")
            {
                hSPf0->Fill(silPoint0.Y(), silPoint0.Z());
            }
            if(layer0 == "l0")
            {
                hSPl0->Fill(silPoint0.X(), silPoint0.Z());
            }
            if(layer0 == "r0")
            {
                hSPr0->Fill(silPoint0.X(), silPoint0.Z());
            }
            hTheta3CM->Fill(theta3CMBefore);        // only thetaCm that enter our cuts
            hTheta3Lab->Fill(theta3Lab * TMath::RadToDeg());
            // write to TTree
            Eex_tree = ExRec;
            theta3CM_tree = theta3CM * TMath::RadToDeg();
            EVertex_tree = T3Rec;
            theta3Lab_tree = theta3Lab * TMath::RadToDeg();
            phi3CM_tree = phi3CM;
            outTree->Fill();
        }
    }

    outFile->Write();
    outFile->Close();

    // Compute efficiency
    auto* effCM {new TEfficiency {*hTheta3CM, *hThetaCMAll}};
    effCM->SetNameTitle("effCM", " #epsilon_{TOT} (#theta_{CM});#epsilon;#theta_{CM} [#circ]");
    auto* effLab {new TEfficiency {*hTheta3Lab, *hThetaLabAll}};
    effLab->SetNameTitle("effLab", "#epsilon_{TOT} (#theta_{Lab});#epsilon;#theta_{Lab} [#circ]");
    // Draw if not running for multiple Exs
    if(inspect)
    {
        auto* cSP {new TCanvas {"cSP", "Sil Points"}};
        cSP->DivideSquare(3);
        cSP->cd(1);
        hSPf0->DrawClone("colz");
        sils->GetLayer("f0").GetSilMatrix()->Draw();
        cSP->cd(2);
        hSPl0->DrawClone("colz");
        sils->GetLayer("l0").GetSilMatrix()->Draw();
        cSP->cd(3);
        hSPr0->DrawClone("colz");
        sils->GetLayer("r0").GetSilMatrix()->Draw();

        auto* c0 {new TCanvas {"c0", "Sim inspect 0"}};
        c0->DivideSquare(6);
        c0->cd(1);
        hKin->DrawClone("colz");
        // Draw theo kin
        kinTheo->SetBeamEnergyAndEx(Tbeam, Ex);
        auto* gtheo {kinTheo->GetKinematicLine3()};
        gtheo->Draw("same");
        c0->cd(2);
        hKinRec->DrawClone("colz");
        // gtheo->Draw("l");
        c0->cd(3);
        hRP_X->DrawClone();
        c0->cd(4);
        hEx->DrawClone();
        c0->cd(5);
        hThetaCMAll->SetTitle("All #theta_{CM}");
        hThetaCMAll->DrawClone();
        c0->cd(6);
        hTheta3CM->SetTitle("#theta_{CM} in cuts");
        hTheta3CM->DrawClone();

        auto* c1 {new TCanvas {"c1", "Sim inspect 1"}};
        c1->DivideSquare(6);
        c1->cd(1);
        hPhi3CM->DrawClone();
        c1->cd(2);
        hPhiAll->DrawClone();

        auto* cEff {new TCanvas {"cEff", "Eff plots"}};
        cEff->DivideSquare(7);
        cEff->cd(2);
        effCM->Draw("apl");
        cEff->cd(3);
        effLab->Draw("apl");
        cEff->cd(4);
        hThetaCMAll->DrawClone();
        cEff->cd(5);
        hTheta3CM->DrawClone();
        cEff->cd(6);
        hThetaLabAll->DrawClone();
        cEff->cd(7);
        hTheta3Lab->DrawClone();
    }
}
#endif
