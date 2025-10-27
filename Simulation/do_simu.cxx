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

#include "./Histos.h"

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
    auto niter {static_cast<int>(1e5)};
    gRandom->SetSeed(0);
    // Initialize detectors
    // TPC
    ActRoot::TPCParameters tpc {"Actar"};
    std::cout << "TPC: " << tpc.X() << " " << tpc.Y() << " " << tpc.Z() << '\n';
    // Silicons
    auto* sils {new ActPhysics::SilSpecs};
    std::string silConfig("silspecs_simu"); // no front silicons, only lateral ones
    sils->ReadFile("../configs/" + silConfig + ".conf");
    sils->Print();
    const double sigmaSil {0.060 / 2.355}; // Si resolution
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
            layer.MoveZTo(75, {3});
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
    bool RestOfBeamLine {false}; // If true enables CFA and mylar of entrance
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
    // srim->ReadTable("beamMylar", path + beam + "_" + Mylar + ".txt");
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
            TbeamRand = srim->SlowWithStraggling("beamMylar", TbeamRand, 0.0038); // Mylar CFA
            TbeamRand = srim->SlowWithStraggling("beamCFA", TbeamRand, 19);       // Gas CFA
            TbeamRand = srim->SlowWithStraggling("beamMylar", TbeamRand, 0.012); // Entrance window of ACTAR
            TbeamRand = srim->SlowWithStraggling("beam", TbeamRand, 48); // Gas before pad plane
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
            if(isThereXS)
            {
                while(theta3CMBefore < 0)
                {
                    // theta3CMBefore = xs->SampleCDF(gRandom->Uniform());
                    theta3CMBefore = xs->SampleHist();
                    // std::cout << theta3CMBefore << std::endl;
                } // sample in deg
            }
            else
            {
                theta3CMBefore = TMath::ACos(gRandom->Uniform(-1, 1)) * TMath::RadToDeg();
            }
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

        // Fill Beam ELoss
        if(vertex.X() > 20)
        {
            double DeltaE {Tbeam - srim->SlowWithStraggling("beam", Tbeam, 20)};
            double E {TbeamRand - TbeamCorr};
            hBeamELoss_E->Fill(E, DeltaE);
        }

        // Fill thetaCMall and kin
        hKin->Fill(theta3LabSampled * TMath::RadToDeg(), T3Lab);
        hThetaCMAll->Fill(theta3CMBefore);
        // Extract direction
        hThetaLabAll->Fill(theta3Lab * TMath::RadToDeg());
        hPhiLab->Fill(phi3Lab * TMath::RadToDeg());
        hPhiCM->Fill(phi3CM * TMath::RadToDeg());
        XYZVector direction {TMath::Cos(theta3Lab), TMath::Sin(theta3Lab) * TMath::Sin(phi3Lab),
                             TMath::Sin(theta3Lab) * TMath::Cos(phi3Lab)};
        // Analysis of heavy particle range, 12Li decays to 11Li

        hTheta4Lab->Fill(theta4Lab * TMath::RadToDeg());

        XYZVector directionHeavy {TMath::Cos(theta4Lab), TMath::Sin(theta4Lab) * TMath::Sin(phi4Lab),
                                  TMath::Sin(theta4Lab) * TMath::Cos(phi4Lab)};
        hT4Lab->Fill(T4Lab);
        double T4LabAfterGas {
            srim->Slow("heavy", T4Lab, (256 - vertex.X()) + 100)}; // Very low angle, so not took into acount
        hT4LabAfterGas->Fill(T4LabAfterGas);
        double rangeHeavyInGas {srim->EvalRange("heavy", T4Lab)};
        hRangeHeavyInGas->Fill(rangeHeavyInGas);
        double rangeHeavyInSilAfterGas {srim->EvalRange("heavyInSil", T4LabAfterGas)};
        hRangeHeavyInSilAfterGas->Fill(rangeHeavyInSilAfterGas);
        // Analysis of unreacted beam particle energy in silicons
        double TbeamAfterGas {srim->Slow("beam", TbeamCorr, 256 + 100)}; // ACTAR lenght + silicon separation
        hTbeamAfterGas->Fill(TbeamAfterGas);
        // Threshold L1, particles that stop in actar. Check before doing the continues
        double rangeInGas {srim->EvalRange("light", T3Lab)};
        ROOT::Math::XYZPoint finalPointGas {vertex + rangeInGas * direction.Unit()};
        if(0 <= finalPointGas.X() && finalPointGas.X() <= 256 && 0 <= finalPointGas.Y() && finalPointGas.Y() <= 256 &&
           0 <= finalPointGas.Z() && finalPointGas.Z() <= 256)
        {
            double distanceXY {
                TMath::Sqrt(pow(vertex.X() - finalPointGas.X(), 2) + pow(vertex.Y() - finalPointGas.Y(), 2))};
            if(distanceXY >= 0)
            {
                hThetaVertexInGas->Fill(theta3Lab * TMath::RadToDeg(), vertex.X());
                hRangeGas->Fill(rangeInGas);
                hKinStopInside->Fill(theta3Lab * TMath::RadToDeg(), T3Lab);

                hRangeInGasSin->Fill(theta3Lab * TMath::RadToDeg(), rangeInGas * TMath::Sin(theta3Lab));
                hRange1->Fill(theta3Lab * TMath::RadToDeg(), rangeInGas);
                hT3test->Fill(T3Lab);

                // Check if the event made out the cuts of an L1 type event
                // CheckL1Acceptance(direction, vertex, finalPointGas, rangeInGas, 8, 16);

                // Fill also eff hists
                if(!onlyLatSils)
                {
                    hThetaCM->Fill(theta3CMBefore);
                    hThetaLab->Fill(theta3Lab * TMath::RadToDeg());
                    hThetaLabStopInside->Fill(theta3Lab * TMath::RadToDeg());
                    hThetaCMStopInside->Fill(theta3CMBefore);
                }
            }
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
        // Check if corresponds to a hit when the detector was off or on
        if(!AcceptHit(silEfficiencies, layer0, silIndex0))
        {
            continue; // if not accepted, go to next iteration
        }
        // get only the events that impact in l0 and r0 for experiment related simulations
        if(onlyLatSils && (layer0 != "l0" && layer0 != "r0"))
            continue;

        // Check now hit of heavy particle
        int silIndexHeavy0 = -1;
        ROOT::Math::XYZPoint silPointHeavy0;
        std::string layerHeavy0;
        std::tie(silIndexHeavy0, silPointHeavy0) = sils->FindSPInLayer("f2", vertex, directionHeavy);
        bool hitHeavy0 {false};
        if(silIndexHeavy0 != -1)
        {
            hitHeavy0 = true;
            layerHeavy0 = "f2";
        }
        // And now check hit in f3
        int silIndexHeavy1 = -1;
        ROOT::Math::XYZPoint silPointHeavy1;
        std::string layerHeavy1;
        bool hitHeavy1 {false};
        if(hitHeavy0)
        {
            std::tie(silIndexHeavy1, silPointHeavy1) = sils->FindSPInLayer("f3", vertex, directionHeavy);
            if(silIndexHeavy1 != -1)
            {
                hitHeavy1 = true;
                layerHeavy1 = "f3";
            }
        }
        // Compute eLoss in silicons for Heavy particle
        if(hitHeavy0 && hitHeavy1)
        {
            auto normalHeavy {sils->GetLayer(layerHeavy0).GetNormal()};
            auto angleWithNormalHeavy {TMath::ACos(directionHeavy.Unit().Dot(normalHeavy.Unit()))};
            auto T4InSil0 {srim->SlowWithStraggling("heavy", T4Lab, (silPointHeavy0 - vertex).R())};
            auto T4AfterSil0 {srim->SlowWithStraggling(
                "heavyInSil", T4InSil0, sils->GetLayer(layerHeavy0).GetUnit().GetThickness(), angleWithNormalHeavy)};
            auto T4InSil1 {srim->SlowWithStraggling("heavy", T4AfterSil0, (silPointHeavy1 - silPointHeavy0).R())};
            auto T4AfterSil1 {srim->SlowWithStraggling(
                "heavyInSil", T4InSil1, sils->GetLayer(layerHeavy1).GetUnit().GetThickness(), angleWithNormalHeavy)};
            auto eLossf2 = T4InSil0 - T4AfterSil0;
            auto eLossf3 = T4InSil1 - T4AfterSil1;
            hHeavyPIDTelescope->Fill(eLossf2, eLossf3);
            // Asign eLoss to ttree data
            eLossSilf2Heavy = eLossf2;
            eLossSilf3Heavy = eLossf3;
        }

        // Fill the SilData for silicon hit on all cases
        // FillSiliconHitsNoCuts(silData, theta3Lab, phi3Lab, T3Lab, theta4Lab, phi4Lab, T4Lab, vertex, AllsilLayers,
        // sils, *silRes, srim); Fill tree with no cuts
        theta3CM_treeNoCut = theta3CMBefore;
        theta3Lab_treeNoCut = theta3Lab * TMath::RadToDeg();
        phi3CM_treeNoCut = phi3CM;
        T3Lab_treeNoCut = T3Lab;
        theta4Lab_treeNoCut = theta4Lab;
        T4Lab_treeNoCut = T4Lab;
        weight_treeNoCut = weight;
        phi4CM_treeNoCut = phi4Lab;
        RPNoCuts_tree = vertex;
        silData_tree = silData;
        outTreeNoCut->Fill();
        // "f0": key name of layer to check for SP
        // silIndex == -1 if NO SP
        // else, returns the silicon index
        // silPoint0: SP in millimetres. in standard ACTAR frame (same as analysis)
        if(silIndex0 == -1)
        {
            hKinNoReachSil->Fill(theta3Lab * TMath::RadToDeg(), T3Lab);
            hThetaLabNoReachSil->Fill(theta3Lab * TMath::RadToDeg());
            hThetaCMNoReachSil->Fill(theta3CMBefore);
            continue;
        }
        hThetaLabDirectionSil->Fill(theta3Lab * TMath::RadToDeg());
        hThetaCMDirectionSil->Fill(theta3CMBefore);

        // Fill hist to get counts in laterals with no cuts
        if(layer0 == "l0")
        {
            hSPl0noCut->Fill(silPoint0.X(), silPoint0.Z());
            hKinLatSil->Fill(theta3Lab * TMath::RadToDeg(), T3Lab);
            hdistanceSPtol0->Fill((silPoint0 - vertex).R());
        }
        if(layer0 == "r0")
        {
            hSPr0noCut->Fill(silPoint0.X(), silPoint0.Z());
            hKinLatSil->Fill(theta3Lab * TMath::RadToDeg(), T3Lab);
            hdistanceSPtor0->Fill((silPoint0 - vertex).R());
        }

        // Slow down light in gas
        auto T3AtSil {srim->SlowWithStraggling("light", T3Lab, (silPoint0 - vertex).R())};
        // Check if stopped
        ApplyNaN(T3AtSil);
        if(std::isnan(T3AtSil))
        {
            hThetaLabNoReachSil->Fill(theta3Lab * TMath::RadToDeg());
            hThetaCMNoReachSil->Fill(theta3CMBefore);
            hRPxEStoppedGas->Fill(vertex.X(), T3Lab);
            if(layer0 == "l0" || layer0 == "r0")
            {
                hTheta3LabNoReachLatSil->Fill(theta3Lab * TMath::RadToDeg());
                hT3LabNoReachLatSil->Fill(T3Lab);
                if(0 <= finalPointGas.X() && finalPointGas.X() <= 256 && 0 <= finalPointGas.Y() &&
                   finalPointGas.Y() <= 256 && 0 <= finalPointGas.Z() && finalPointGas.Z() <= 256)
                {
                    hKinLatStopInside->Fill(theta3Lab * TMath::RadToDeg(), T3Lab);
                }
                else
                {
                    hKinStopBetweenSilAndGas->Fill(theta3Lab * TMath::RadToDeg(), T3Lab);
                }
            }
            continue;
        }
        double DeltaELength {-1};
        if(layer0 == "f0")
        {
            auto line {
                new ActRoot::Line(ActRoot::CastXYZPoint<float>(vertex), ActRoot::CastXYZVector<float>(direction), 0)};
            auto boundaryPoint {line->MoveToX(256)};
            auto T3AtSilLight {T3Lab - srim->SlowWithStraggling("light", T3Lab, (boundaryPoint - vertex).R())};
            DeltaELength = T3AtSilLight / (boundaryPoint - vertex).R() * 1000;
            hDeltaEGasThetaLight->Fill(theta3Lab * TMath::RadToDeg(), DeltaELength);
            // if(DeltaELength > 2)
            //     hKinDebug->Fill(theta3Lab * TMath::RadToDeg(), T3Lab);
            delete line;
        }
        if(layer0 == "l0")
        {
            auto line {
                new ActRoot::Line(ActRoot::CastXYZPoint<float>(vertex), ActRoot::CastXYZVector<float>(direction), 0)};
            auto boundaryPoint {line->MoveToY(256)};
            auto T3AtSilLight {T3Lab - srim->SlowWithStraggling("light", T3Lab, (boundaryPoint - vertex).R())};
            DeltaELength = T3AtSilLight / (boundaryPoint - vertex).R() * 1000;
            hDeltaEGasThetaLightSide->Fill(theta3Lab * TMath::RadToDeg(), DeltaELength);
            delete line;
            // if(DeltaELength > 2)
            //     hKinDebug->Fill(theta3Lab * TMath::RadToDeg(), T3Lab);
        }
        if(layer0 == "r0")
        {
            auto line {
                new ActRoot::Line(ActRoot::CastXYZPoint<float>(vertex), ActRoot::CastXYZVector<float>(direction), 0)};
            auto boundaryPoint {line->MoveToY(0)};
            auto T3AtSilLight {T3Lab - srim->SlowWithStraggling("light", T3Lab, (boundaryPoint - vertex).R())};
            DeltaELength = T3AtSilLight / (boundaryPoint - vertex).R() * 1000;
            hDeltaEGasThetaLightSide->Fill(theta3Lab * TMath::RadToDeg(), DeltaELength);
            delete line;
            // if(DeltaELength > 2)
            //     hKinDebug->Fill(theta3Lab * TMath::RadToDeg(), T3Lab);
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
            hThetaLabNoReachSil->Fill(theta3Lab * TMath::RadToDeg());
            hThetaCMNoReachSil->Fill(theta3CMBefore);

            hTheta3LabNoReachLatSil->Fill(theta3Lab * TMath::RadToDeg());
            hT3LabNoReachLatSil->Fill(T3Lab);

            if(layer0 == "l0" || layer0 == "r0")
                hKinStopBetweenSilAndGas->Fill(
                    theta3Lab * TMath::RadToDeg(),
                    T3Lab); // No having more energy than threshold = stop between gas and silicon
            continue;
        }

        // Analysis of punchthrough in lateral silicons
        if(T3AfterSil0 > 0 && (layer0 == "l0" || layer0 == "r0"))
        {
            hTheta3LabLatSilPunch->Fill(theta3Lab * TMath::RadToDeg());
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
                hRPxEfirstSil->Fill(vertex.X(), T3Rec);
            }
            else
            {
                auto T3Rec0 {srim->EvalInitialEnergy("light", eLoss0, (silPoint0 - vertex).R())};
                hRPxEfirstSil->Fill(vertex.X(),
                                    T3Rec0); // This is to show what would be obtained if there was just one sil

                // Reconstruction ob T3 with 2 silicon layers
                auto T3Rec1 {srim->EvalInitialEnergy("light", eLoss1, (silPoint1 - silPoint0).R())};
                T3Rec = srim->EvalInitialEnergy("light", eLoss0 + T3Rec1, (silPoint0 - vertex).R());
            }

            auto ExRec {kin->ReconstructExcitationEnergy(T3Rec, theta3Lab)};

            // Fill
            hKinRec->Fill(theta3Lab * TMath::RadToDeg(), T3Rec); // after reconstruction
            hEx->Fill(ExRec, weight * alpha);
            hRP->Fill(vertex.X(), vertex.Y());
            hRP_X->Fill(vertex.X());
            hRP_ZY->Fill(vertex.Y(), vertex.Z());
            if(layer0 == "f0")
            {
                hSPf0->Fill(silPoint0.Y(), silPoint0.Z());
            }
            if(layer0 == "l0")
            {
                hSPl0->Fill(silPoint0.X(), silPoint0.Z());
                hTheta3LabReachLatSil->Fill(theta3Lab * TMath::RadToDeg());
                hKinLatSilReach->Fill(theta3Lab * TMath::RadToDeg(), T3Lab);
                hEVertexESil->Fill(T3Lab, T3AtSil);
            }
            if(layer0 == "r0")
            {
                hSPr0->Fill(silPoint0.X(), silPoint0.Z());
                hTheta3LabReachLatSil->Fill(theta3Lab * TMath::RadToDeg());
                hKinLatSilReach->Fill(theta3Lab * TMath::RadToDeg(), T3Lab);
                hEVertexESil->Fill(T3Lab, T3AtSil);
            }
            hThetaCMSilicon->Fill(theta3CMBefore); // Other hThetaCM has efficiency for L1 events also
            hThetaCM->Fill(theta3CMBefore);        // only thetaCm that enter our cuts
            hThetaLab->Fill(theta3Lab * TMath::RadToDeg());
            hThetaCMThetaLab->Fill(theta3CMBefore, theta3LabSampled * TMath::RadToDeg());
            hRPxEbothSil->Fill(vertex.X(), T3Rec);
            hRPeffIn->Fill(vertex.X());

            hThetaLabMeassureSil->Fill(theta3Lab * TMath::RadToDeg());
            hThetaCMMeassureSil->Fill(theta3CMBefore);

            // Debufg Plots
            hKinDebug->Fill(theta3Lab * TMath::RadToDeg(), T3Rec);
            hKinHeavyDebug->Fill(theta4Lab * TMath::RadToDeg(), T4Lab);
            hTheta3Theta4LabDebug->Fill(theta3Lab * TMath::RadToDeg(), theta4Lab * TMath::RadToDeg());
            if(!(silIndexHeavy1 == -1))
            {
                hPIDHeavyDebug->Fill(eLossSilf2Heavy, eLossSilf3Heavy);
            }
            hSPf0HeavyDebug->Fill(silPointHeavy0.Y(), silPointHeavy0.Z());

            // write to TTree
            Eex_tree = ExRec;
            theta3CM_tree = theta3CM * TMath::RadToDeg();
            EVertex_tree = T3Rec;
            theta3Lab_tree = theta3Lab * TMath::RadToDeg();
            phi3CM_tree = phi3CM;
            DeltaEgas_tree = DeltaELength;
            DeltaESil_tree = eLoss0;
            DeltaESil1_tree = eLoss1;
            outTree->Fill();
            theta4Lab_tree = theta4Lab;
            phi4Lab_tree = phi4Lab;
            T4Lab_tree = T4Lab;
            RP_tree = vertex;
            weight_tree = weight;
            outTreeHeavy->Fill();
        }
        delete silData; // delete silData to avoid memory leaks
    }

    outFile->Write();
    outFile->Close();

    // Compute efficiency
    auto* effCM {new TEfficiency {*hThetaCM, *hThetaCMAll}};
    effCM->SetNameTitle("effCM", " #epsilon_{TOT} (#theta_{CM});#epsilon;#theta_{CM} [#circ]");
    auto* effLab {new TEfficiency {*hThetaLab, *hThetaLabAll}};
    effLab->SetNameTitle("effLab", "#epsilon_{TOT} (#theta_{Lab});#epsilon;#theta_{Lab} [#circ]");
    auto* effLatSil {new TEfficiency {*hTheta3LabReachLatSil, *hThetaLabAll}};
    effLatSil->SetNameTitle("effLatSil", "#epsilon_{TOT} (#theta_{Lab}) Lateral Silicon;#epsilon;#theta_{Lab} [#circ]");
    auto* effStopInside {new TEfficiency {*hThetaLabStopInside, *hThetaLabAll}};
    effStopInside->SetNameTitle("effStopInside",
                                "#epsilon_{TOT} (#theta_{Lab}) events stopped inside;#epsilon;#theta_{Lab} [#circ]");
    auto effDetectionSilLab {new TEfficiency {*hThetaLabMeassureSil, *hThetaLabDirectionSil}};
    effDetectionSilLab->SetNameTitle("effDetectionSil", "#epsilon_{Sil} (#theta_{Lab}) events detected in silicons "
                                                        "that go to them;#epsilon;#theta_{Lab} [#circ]");
    auto effDetectionSilCM {new TEfficiency {*hThetaCMMeassureSil, *hThetaCMDirectionSil}};
    effDetectionSilCM->SetNameTitle("effDetectionSilCM", "#epsilon_{Sil} (#theta_{CM}) events detected in silicons "
                                                         "that go to them;#epsilon;#theta_{CM} [#circ]");
    auto effDetectionL1Lab {new TEfficiency {*hThetaLabStopInside, *hThetaLabNoReachSil}};
    effDetectionL1Lab->SetNameTitle("effDetectionL1Lab", "#epsilon_{L1} (#theta_{Lab}) events detected L1 (no Sil "
                                                         "signal);#epsilon;#theta_{Lab} [#circ]");
    auto effDetectionL1CM {new TEfficiency {*hThetaCMStopInside, *hThetaCMNoReachSil}};
    effDetectionL1CM->SetNameTitle("effDetectionL1CM", "#epsilon_{L1} (#theta_{CM}) events detected L1 (no Sil "
                                                       "signal);#epsilon;#theta_{CM} [#circ]");
    // Now efficiency in intervals
    auto* effIntervals {new TEfficiency {*hRPeffIn, *hRPeffAll}};
    effIntervals->SetNameTitle("effIntervals", "#RPx efficiency;#epsilon;#RPx [#mm]");

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
        hRP->DrawClone("colz");
        c0->cd(5);
        hThetaCMAll->SetTitle("All #theta_{CM}");
        hThetaCMAll->DrawClone();
        c0->cd(6);
        hThetaCM->SetTitle("#theta_{CM} in cuts");
        hThetaCM->DrawClone();

        auto* c1 {new TCanvas {"c1", "Sim inspect 1"}};
        c1->DivideSquare(6);
        c1->cd(1);
        hPhiLab->DrawClone();
        c1->cd(2);
        hEx->DrawClone();
        c1->cd(3);
        hThetaCMThetaLab->DrawClone("colz");
        // auto* gCMLab {kin->GetThetaLabvsThetaCMLine()};
        // gCMLab->Draw("l");
        c1->cd(4);
        hRPxEfirstSil->DrawClone();
        c1->cd(5);
        hRPxEbothSil->DrawClone();
        c1->cd(6);
        hRPxEStoppedGas->DrawClone();

        auto* cEff {new TCanvas {"cEff", "Eff plots"}};
        cEff->DivideSquare(7);
        cEff->cd(1);
        effIntervals->Draw("apl");
        cEff->cd(2);
        effCM->Draw("apl");
        cEff->cd(3);
        effLab->Draw("apl");
        cEff->cd(4);
        hThetaCMAll->DrawClone();
        cEff->cd(5);
        hThetaCM->DrawClone();
        cEff->cd(6);
        hThetaLabAll->DrawClone();
        cEff->cd(7);
        hThetaLab->DrawClone();

        auto* cL1 {new TCanvas {"cL1", "Plots for L1 events"}};
        cL1->DivideSquare(6);
        cL1->cd(1);
        hThetaLabStopInside->DrawClone();
        cL1->cd(2);
        hThetaVertexInGas->DrawClone("colz");
        cL1->cd(3);
        hRangeGas->DrawClone();

        auto* cHeavy {new TCanvas {"cHeavy", "Heavy Particle Study"}};
        cHeavy->DivideSquare(6);
        cHeavy->cd(1);
        hT4Lab->DrawClone();
        cHeavy->cd(2);
        hRangeHeavyInGas->DrawClone();
        cHeavy->cd(3);
        hRangeHeavyInSilAfterGas->DrawClone();
        cHeavy->cd(4);
        hTheta4Lab->DrawClone();
        cHeavy->cd(5);
        hT4LabAfterGas->DrawClone();
        cHeavy->cd(6);
        hTbeamAfterGas->DrawClone();

        auto* cLat {new TCanvas {"cLatnoCut", "lateral Silicon Study"}};
        cLat->DivideSquare(7);
        cLat->cd(1);
        hSPl0noCut->DrawClone("colz");
        sils->GetLayer("l0").GetSilMatrix()->DrawClone();
        cLat->cd(2);
        hSPr0noCut->DrawClone("colz");
        sils->GetLayer("r0").GetSilMatrix()->DrawClone();
        cLat->cd(3);
        hTheta3LabReachLatSil->DrawClone();
        cLat->cd(4);
        hTheta3LabNoReachLatSil->DrawClone();
        cLat->cd(5);
        hT3LabNoReachLatSil->DrawClone();
        cLat->cd(6);
        hTheta3LabLatSilPunch->DrawClone();
        cLat->cd(7);
        hEVertexESil->DrawClone("colz");
        cLat->cd(8);
        hdistanceSPtol0->DrawClone();

        auto cLatKin {new TCanvas {"cLatKin", "Kinematics for lateral silicons"}};
        cLatKin->DivideSquare(10);
        cLatKin->cd(1);
        hKinLatSil->DrawClone();
        cLatKin->cd(2);
        hKinLatSilReach->DrawClone();
        cLatKin->cd(3);
        hKinStopBetweenSilAndGas->DrawClone();
        cLatKin->cd(4);
        hKinLatStopInside->DrawClone();
        cLatKin->cd(5);
        hKinStopInside->DrawClone();
        cLatKin->cd(6);
        hKinNoReachSil->DrawClone();
        cLatKin->cd(7);
        hTheta3LabReachLatSil->DrawClone();
        cLatKin->cd(8);
        hThetaLabStopInside->DrawClone();
        cLatKin->cd(9);
        effLatSil->Draw("apl");
        cLatKin->cd(10);
        effStopInside->Draw("apl");

        auto cProtonStopped {new TCanvas {"cProtonStopped", "Protons stopped in gas"}};
        cProtonStopped->DivideSquare(4);
        cProtonStopped->cd(1);
        hRangeInGasSin->DrawClone();
        // hRange1->DrawClone("same");
        cProtonStopped->cd(2);
        hKinStopInside->DrawClone();
        // gtheo->Draw("l");
        cProtonStopped->cd(3);
        hT3test->DrawClone();

        auto cCheckAsymmetry {new TCanvas {"cCheckAsymmetry", "Check Assymetry"}};
        cCheckAsymmetry->DivideSquare(6);
        cCheckAsymmetry->cd(1);
        hPhiCM->DrawClone();
        cCheckAsymmetry->cd(2);
        hPhiLab->DrawClone();
        cCheckAsymmetry->cd(3);
        hRP->DrawClone("colz");
        cCheckAsymmetry->cd(4);
        hRP_ZY->DrawClone("colz");
        cCheckAsymmetry->cd(5);
        hdistanceSPtol0->DrawClone();
        cCheckAsymmetry->cd(6);
        hdistanceSPtor0->DrawClone();

        auto cEfficiency {new TCanvas {"cEfficiency", "Efficiency"}};
        cEfficiency->DivideSquare(6);
        cEfficiency->cd(1);
        effLab->Draw("apl");
        cEfficiency->cd(2);
        effDetectionSilLab->Draw("apl");
        cEfficiency->cd(3);
        effDetectionL1Lab->Draw("apl");
        cEfficiency->cd(4);
        effCM->Draw("apl");
        cEfficiency->cd(5);
        effDetectionSilCM->Draw("apl");
        cEfficiency->cd(6);
        effDetectionL1CM->Draw("apl");

        auto cBeam {new TCanvas {"cBeam", "Beam Spectra"}};
        cBeam->DivideSquare(4);
        cBeam->cd(1);
        hTbeam->DrawClone();

        auto cDebug {new TCanvas {"cDebug", "Debugging plots"}};
        cDebug->DivideSquare(6);
        cDebug->cd(1);
        hKinDebug->DrawClone("colz");
        cDebug->cd(2);
        hKinHeavyDebug->DrawClone("colz");
        cDebug->cd(3);
        hTheta3Theta4LabDebug->DrawClone("colz");
        // cDebug->cd(4);
        // hPIDHeavyDebug->DrawClone("colz");
        cDebug->cd(5);
        hSPf0HeavyDebug->DrawClone("colz");

        auto cBeamCharge = new TCanvas("cBeamCharge", "Beam energy loss");
        hBeamELoss_E->DrawClone("colz");
    }
}
#endif
