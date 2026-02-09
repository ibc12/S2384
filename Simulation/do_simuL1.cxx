#ifndef triumf_all_cxx
#define triumf_all_cxx
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
#include <random>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>

#include "../PostAnalysis/HistConfig.h"

using XYZPoint = ROOT::Math::XYZPoint;
using XYZVector = ROOT::Math::XYZVector;

struct BeamOffset
{
    double offset;   // mm
    double fraction; // must sum to 1
};

using padPlane = std::map<std::pair<int, int>, double>;

using BeamOffsetMap = std::map<std::string, std::vector<BeamOffset>>;

std::pair<XYZPoint, XYZPoint> SampleVertex(double meanZ, double sigmaZ, TH3D* h, double lengthX)
{

    // X is always common for both manners
    double Xstart {0};
    double Xrp {gRandom->Uniform() * lengthX};
    // Y depends completely on the method of calculation
    double Ystart {-1};
    double Yrp {-1};
    // Z of beam at entrance
    double Zstart {gRandom->Gaus(meanZ, sigmaZ)};
    double Zrp {-1};
    // Ystart in this case is sampled from the histogram itself!
    double thetaXY {};
    double thetaXZ {};
    h->GetRandom3(Ystart, thetaXY, thetaXZ);
    // Mind that Y is not centred in the histogram value!
    // Rp values are computed as follows:
    Yrp = Ystart - Xrp * TMath::Tan(thetaXY * TMath::DegToRad());
    Zrp = Zstart - Xrp * TMath::Tan(thetaXZ * TMath::DegToRad());
    XYZPoint start {Xstart, Ystart, Zstart};
    XYZPoint vertex {Xrp, Yrp, Zrp};
    return {std::move(start), std::move(vertex)};
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
    double sigma {1.5 / 2.355}; // FWHM to sigma
    theta = gRandom->Gaus(theta, sigma * TMath::DegToRad());
}

double RandomizeBeamEnergy(double Tini, double sigma)
{
    return gRandom->Gaus(Tini, sigma);
}

// ============================================================
// Polya function for gain distribution in Micromegas
// ============================================================
double Polya(double Gmean, double theta)
{
    static std::mt19937 gen(std::random_device {}());
    std::gamma_distribution<double> gamma(theta + 1.0, Gmean / (theta + 1.0));
    return gamma(gen);
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
        else if(Ex == 0.477)
        {
            isThereXS = true;
            TString data_to_read {TString::Format("./Inputs/xs/%s/dd/g1.dat", beam.c_str())};
            xs->ReadFile(data_to_read.Data());
            std::cout << "Total xs: " << xs->GetTotalXSmbarn() << std::endl;
        }
    }
    return isThereXS;
}

//     // ============================================================
//     // Add charge to voxel map
//     // ============================================================
//     void AddChargeToVoxel(double x, double y, double z, double charge, std::map<voxelKey, ActRoot::Voxel>& voxelMap)
//     {
//         int ix = std::floor(x / voxelSize);
//         int iy = std::floor(y / voxelSize);
//         int iz = std::floor(z / voxelSize);
//     
//         voxelKey key {ix, iy, iz};
//     
//         auto it = voxelMap.find(key);
//         if(it == voxelMap.end())
//         {
//             ActRoot::Voxel v;
//             v.SetPosition(ActRoot::Voxel::XYZPointF {float((ix + 0.5) * voxelSize), float((iy + 0.5) * voxelSize),
//                                                      float((iz + 0.5) * voxelSize)});
//             v.SetCharge(charge);
//             voxelMap.emplace(key, v);
//         }
//         else
//         {
//             it->second.SetCharge(it->second.GetCharge() + charge);
//         }
//     }
//     
//     // ============================================================
//     // Divide segment → electrons → voxels
//     // ============================================================
//     void DivideSegmentInPortions(double eLoss, int nPortions, const XYZPoint& center,
//                                  std::map<voxelKey, ActRoot::Voxel>& voxelMap, std::vector<XYZPoint>& electrons)
//     {
//         if(eLoss <= 0 || nPortions <= 0)
//             return;
//     
//         const double W = 36.4 * 0.95 + 34.3 * 0.05; // eV / electron in 95% D2 + 5% CF4
//         const double diffT = 0.10;                  // mm / sqrt(mm)
//         const double diffL = 0;                     // mm / sqrt(mm) No diffusion in z direction
//     
//         double h = center.Y();
//         if(h <= 0)
//             return;
//     
//         double sigmaT = diffT * std::sqrt(h);
//         double sigmaL = diffL * std::sqrt(h);
//     
//         double portionE = eLoss / nPortions; // MeV
//         double meanNe = portionE * 1e6 / W;
//     
//         for(int i = 0; i < nPortions; i++)
//         {
//             int nElectrons = gRandom->Poisson(meanNe);
//             if(nElectrons <= 0)
//                 continue;
//     
//             for(int e = 0; e < nElectrons; e++)
//             {
//                 double x = gRandom->Gaus(center.X(), sigmaT);
//                 double y = gRandom->Gaus(center.Y(), sigmaT);
//                 double z = gRandom->Gaus(center.Z(), sigmaL);
//     
//                 double gain = Polya(Gmean, theta);
//                 AddChargeToVoxel(x, y, z, gain, voxelMap);
//                 electrons.emplace_back(x, y, z);
//             }
//         }
//     }
//     
//     // ============================================================
//     // Divide track using SRIM
//     // ============================================================
//     void DivideTrackInSegments(ActPhysics::SRIM* srim, double range, const XYZVector& dirIn, const XYZPoint& rp,
//                                double step, int nSub, std::map<voxelKey, ActRoot::Voxel>& voxelMap,
//                                std::vector<XYZPoint>& electrons)
//     {
//         XYZVector dir = dirIn.Unit();
//         double E = srim->EvalInverse("light", range);
//     
//         for(double r = 0; r < range; r += step)
//         {
//             double Epost = srim->Slow("light", E, step);
//             double eLoss = E - Epost;
//             E = Epost;
//     
//             if(eLoss <= 0)
//                 continue;
//     
//             XYZPoint center = rp + dir * (r + 0.5 * step);
//             if(center.Y() <= 0)
//                 continue;
//     
//             DivideSegmentInPortions(eLoss, nSub, center, voxelMap, electrons);
//         }
//     }

padPlane
DivideTrackInSegments(ActPhysics::SRIM* srim, double range, XYZVector dir, XYZPoint rp, double step, int nSubSegments)
{
    // Divide track in segments of length = step
    // Compute energy loss in that step
    // Divide that eLoss in the segment in a number of portions

    padPlane padMap {};

    // Normalize the drection
    dir = dir.Unit();
    // Initial energy
    double Eiter {srim->EvalInverse("light", range)};

    for(double r = 0; r < range; r += step)
    {
        double EpostSlow {srim->Slow("light", Eiter, step)};
        double eLoss {Eiter - EpostSlow};

        if(eLoss <= 0)
            continue;

        // Mean height of the segment for the drift
        XYZPoint posSegment = rp + dir * (r + 0.5 * step);
        double hSegment = posSegment.Z();

        // Divide in subsegments to drift
        DivideSegmentInPortions(eLoss, nSubSegments, hSegment, padMap);

        Eiter = EpostSlow;
    }

    return padMap;
}

void DivideSegmentInPortions(double eLoss, int nPortions, double hSegment, padPlane& padMap)
{
    // Each segment divide it in many portions

    // Each portion has the same amount of energy

    // Transform energy in electrons

    // Randomize the electron numbers with Poisson distribution

    // Drift Electrons into pad plane with gausian distribution with width the difusion taking into acount the heigth of
    // the charge

    // Fill map of <pad,pad> , charge that represents the pad plane

    if(eLoss <= 0 || nPortions <= 0)
        return;

    // Constantes
    const double W = 30.0;          // eV/electron
    const double driftDiffT = 0.1;  // mm/sqrt(mm), transversal diffusion
    const double driftDiffL = 0.05; // mm/sqrt(mm), longitudinal diffusion
    const double padSize = 2.0;     // mm

    double portionEnergy = eLoss / nPortions;
    double meanElectrons = portionEnergy * 1e6 / W; // eLoss en MeV a eV a electrones

    for(int i = 0; i < nPortions; i++)
    {
        // Randomize electrons in each portion
        int nElectrons = gRandom->Poisson(meanElectrons);
        if(nElectrons <= 0)
            continue;

        // Difusión
        double sigmaT = driftDiffT * std::sqrt(hSegment); // mm
        double sigmaL = driftDiffL * std::sqrt(hSegment); // mm

        for(int e = 0; e < nElectrons; e++)
        {
            // posición de cada electrón en plano de pads (x,z) con difusión gaussiana
            double xPad = gRandom->Gaus(0.0, sigmaT);
            double zPad = gRandom->Gaus(0.0, sigmaL);

            // Convertir posición a número de pad (0–2 → 0, 2–4 → 1, ...)
            int col = int(xPad / padSize);
            int row = int(zPad / padSize);

            padMap[std::make_pair(row, col)] += 1.0; // acumular carga
        }
    }
    // Check ActSim, to see if the drift is done electron by electron
    // Check If the drift sigma is well computed with the height
}

bool CheckExclusionZone(padPlane padPlane, int nPadThreshold, int yMin, int yMax, int chargeThreshold = 0)
{
    // Check the number of pads outside the exclusion zone defined by yMin and yMax have charge above threshold

    // I have to check, because maybe the consition is not only pads with charge, Thomas explained this in July
}

void do_simu(const std::string& beam, const std::string& target, const std::string& light, const std::string& heavy,
             int neutronPS, int protonPS, double Tbeam, double Ex, bool inspect, int thread = -1)
{
    // set batch mode if inspect is false
    if(!inspect)
        gROOT->SetBatch(true);
    // Set whether is PS or not
    bool isPS {(neutronPS > 0) || (protonPS > 0)};
    // Set number of iterations
    const int niter {static_cast<int>(inspect ? 1e7 : (isPS ? 1e8 : 1e8))};
    gRandom->SetSeed(0);
    // Runner: contains utility functions to execute multiple actions as rotate directions
    ActSim::Runner runner(nullptr, nullptr, gRandom, 0);
    // Initialize detectors
    // TPC
    ActRoot::TPCParameters tpc {"Actar"};
    std::cout << "TPC: " << tpc.X() << " " << tpc.Y() << " " << tpc.Z() << '\n';
    // Vertex sampling and beam z variables
    std::string emittancefilename {};
    std::string beamPositionTimefilename {};
    if(beam == "7Li")
    {
        emittancefilename = {"../Macros/Emittance/Outputs/histos" + beam + ".root"};
        beamPositionTimefilename = ""; // 7Li did not change position during experiment
    }

    else if(beam == "11Li")
    {
        emittancefilename = {"../Macros/Emittance/Outputs/histos" + beam + "_pre.root"};
        beamPositionTimefilename =
            "./Inputs/Efficiencies/beamEmittancePeriods_11Li.dat"; // File with proportion of time that beam was in
                                                                   // each position
    }
    auto beamfile {std::make_unique<TFile>(emittancefilename.c_str())};
    auto* hBeam {beamfile->Get<TH3D>("h3d")};
    if(!hBeam)
        throw std::runtime_error("Could not load beam emittance histogram");
    hBeam->SetDirectory(nullptr);
    beamfile.reset();
    // No silicons for L1

    const double zMeanEntrance {
        135}; // Mean position of ACTAR entrance (exp values z position are not refered to a certain point)
    const double zVertexSigma {0.81}; // From emitance study / always  around the same)

    // No experimental cuts for L1 yet

    // Sigmas
    const double sigmaPercentBeam {0.0019}; // 0,19% beam energy spread (meassured by operators)
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
    // isThereXS = true;
    // xs->ReadFile("../Fits/7Li_dp/Inputs/gs_Daehnik_Delaroche_myself/21.g0");
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
    auto hKin {HistConfig::KinEl.GetHistogram()};
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
    // Reconstructed histos
    auto hEx {HistConfig::Ex.GetHistogram()};
    auto hKinRec {HistConfig::Kin.GetHistogram()};
    hKinRec->SetTitle("Reconstructed Kinetic Energy;#theta_{Lab} [#circ];E_{Vertex} [MeV]");
    auto hRP_X {HistConfig::RPx.GetHistogram()};
    hRP_X->SetTitle("Reconstructed RP;X [mm];Counts");
    auto hRP {HistConfig::RP.GetHistogram()};
    hRP->SetTitle("Reconstructed RP;RP [mm];Counts");

    // Allow multiple theads
    std::string tag {""};
    if(thread > 0)
        tag = "_" + std::to_string(thread);

    // File to save data
    TString fileName {TString::Format("./Outputs/%s/%s_%s_TRIUMF_Eex_%.3f_nPS_%d_pPS_%d%s_L1.root", beam.c_str(),
                                      target.c_str(), light.c_str(), Ex, neutronPS, protonPS, tag.c_str())};
    auto outFile {new TFile(fileName, inspect ? "read" : "recreate")};
    auto* outTree {new TTree("SimulationTTree", "A TTree containing only our Eex obtained by simulation")};
    if(inspect)
        outTree->SetDirectory(nullptr);
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
    double weight_tree {};
    outTree->Branch("weight", &weight_tree);

    // ---- SIMU STARTS HERE ----
    ROOT::EnableImplicitMT();

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
        // if(it >= nextPrint)
        // {
        //     percent = 100 * (it + 1) / niter;
        //     int nchar {percent / percentPrint};
        //     std::cout << "\r" << std::string((int)(percent / percentPrint), '|') << percent << "%";
        //     std::cout.flush();
        //     nextPrint += step;
        // }
        // Sample vertex position
        auto [start, vertex] {SampleVertex(zMeanEntrance, zVertexSigma, hBeam, tpc.X())};
        auto distToVertex {(vertex - start).R()};

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
            TbeamRand = srim->SlowWithStraggling("beamMylar", TbeamRand, 0.0168); // All mylar
            TbeamRand = srim->SlowWithStraggling(
                "beam", TbeamRand, 60); // Gas before pad plane (approximation, not taking into account the angle)
        }
        auto TbeamCorr {srim->SlowWithStraggling("beam", TbeamRand, distToVertex)};
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

        // Propagate track from vertex to silicon wall using SilSpecs class
        // And using the angle with the uncertainty already in
        ROOT::Math::XYZVector dirBeamFrame {TMath::Cos(theta3Lab), TMath::Sin(theta3Lab) * TMath::Sin(phi3Lab),
                                            TMath::Sin(theta3Lab) * TMath::Cos(phi3Lab)};
        ROOT::Math::XYZVector heavyBeamFrame {TMath::Cos(theta4Lab), TMath::Sin(theta4Lab) * TMath::Sin(phi4Lab),
                                              TMath::Sin(theta4Lab) * TMath::Cos(phi4Lab)};
        // Declare beam direction
        auto beamDir {(vertex - start).Unit()};
        // Rotate to world = geometry frame
        auto dirWorldFrame {runner.RotateToWorldFrame(dirBeamFrame, beamDir)};
        auto heavyWorldFrame {runner.RotateToWorldFrame(heavyBeamFrame, beamDir)};


        // Extract direction
        XYZVector direction {TMath::Cos(theta3Lab), TMath::Sin(theta3Lab) * TMath::Sin(phi3Lab),
                             TMath::Sin(theta3Lab) * TMath::Cos(phi3Lab)};
        // L1 condition, particles that stop in actar. Check if track stays inside
        double rangeInGas {srim->EvalRange("light", T3Lab)};
        ROOT::Math::XYZPoint finalPointGas {vertex + rangeInGas * dirWorldFrame.Unit()};
        bool isL1 {0 <= finalPointGas.X() && finalPointGas.X() <= tpc.X() && 0 <= finalPointGas.Y() &&
                   finalPointGas.Y() <= tpc.Y() && 0 <= finalPointGas.Z() && finalPointGas.Z() <= tpc.Z()};
        if(!isL1)
            continue;

        // Exclusion zone from pad 55 to 70 (pads start at 0, so from 108 to 142 mm in Y)

        // Reconstruct track and charge drift, diffusion and deposition


        // Reconstruct Ex!
        bool isOk {};          // no punchthrouhg
        bool cutELoss0 {true}; // for f0 not yet implemented the graphical cuts
        if(isOk && cutELoss0)
        {
            double T3Rec {0};
            double ExRec {0};
            // Fill
            hKinRec->Fill(theta3Lab * TMath::RadToDeg(), T3Rec); // after reconstruction
            hEx->Fill(ExRec, weight);                            // To get real counts weigth * alpha
            hRP_X->Fill(vertex.X());
            hRP->Fill(vertex.X(), vertex.Y());
            hTheta3CM->Fill(theta3CMBefore); // only thetaCm that enter our cuts
            hTheta3Lab->Fill(theta3Lab * TMath::RadToDeg());
            // write to TTree
            Eex_tree = ExRec;
            theta3CM_tree = theta3CM * TMath::RadToDeg();
            EVertex_tree = T3Rec;
            theta3Lab_tree = theta3Lab * TMath::RadToDeg();
            phi3CM_tree = phi3CM;
            weight_tree = weight;
            outTree->Fill();
        }
    }


    // Compute efficiency side, front and total
    auto* effCM {new TEfficiency {*hTheta3CM, *hThetaCMAll}};
    effCM->SetNameTitle("effCM", " #epsilon_{TOT} (#theta_{CM});#epsilon;#theta_{CM} [#circ]");
    auto* effLab {new TEfficiency {*hTheta3Lab, *hThetaLabAll}};
    effLab->SetNameTitle("effLab", "#epsilon_{TOT} (#theta_{Lab});#epsilon;#theta_{Lab} [#circ]");

    // SAVING
    if(!inspect)
    {
        outFile->cd();
        outTree->Write();
        effCM->Write();
        effLab->Write();
        hRP->Write("hRP");
        outFile->Close();
        delete outFile;
        outFile = nullptr;
    }
    // Draw if not running for multiple Exs
    if(inspect)
    {
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
        // Fit hEx
        hEx->Fit("gaus", "", "", Ex - 1, Ex + 1);
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

    // deleting news
    delete srim;
    if(isThereXS)
        delete xs;
}
#endif
