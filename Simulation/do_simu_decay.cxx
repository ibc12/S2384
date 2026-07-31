#ifndef do_simu_decay_cxx
#define do_simu_decay_cxx
#include "ActColors.h"
#include "ActConstants.h"
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
#include "TParameter.h"
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

#include "../PostAnalysis/HistConfig.h"
#include "./do_simu.cxx"

using XYZPoint = ROOT::Math::XYZPoint;
using XYZVector = ROOT::Math::XYZVector;

using BeamOffsetMap = std::map<std::string, std::vector<BeamOffset>>;


// ==========================================================================
// Result of propagating a charge particle (alfa or triton) through one of the two front dE-E telescopes: (f0,f1) or
// (f2,f3).
// ==========================================================================
struct FrontHitResult
{
    bool detected {false};
    int telescopeCode {0}; // 0 = no detected, 1 = f0/f1, 2 = f2/f3
    std::string layerThin {};
    std::string layerThick {};
    int silIndexThin {-1};
    bool punchthrough {false};
    double energy {}; // eLoss total reconstructed (thin [+ thick if punchthrough])
    double angle {};  // thetaLab (in deg)
    ROOT::Math::XYZPoint point {};
};

// Check first the f0/f1 telescope and then the f2/f3; return the first one that gives a valid hit (passes threshold and
// efficiency). Replicate the same pattern of SlowWithStraggling + resolution + ApplyNaN that is already used for
// 'light', generalized to two dE-E pairs instead of just f0->f1.
FrontHitResult TrackFrontWall(const std::string& srimKey, double Tinit, const ROOT::Math::XYZVector& dirWorldFrame,
                              const ROOT::Math::XYZPoint& vertex, ActPhysics::SilSpecs* sils, ActPhysics::SRIM* srim,
                              TF1* silRes, double sigmaFront, std::map<std::string, double>& silEfficiencies)
{
    FrontHitResult res;
    static const std::vector<std::pair<std::string, std::string>> telescopes {{"f0", "f1"}, {"f2", "f3"}};
    for(int t = 0; t < (int)telescopes.size(); t++)
    {
        const auto& [thin, thick] = telescopes[t];
        int silIndexThin {-1};
        ROOT::Math::XYZPoint pointThin;
        std::tie(silIndexThin, pointThin) = sils->FindSPInLayer(thin, vertex, dirWorldFrame);
        if(silIndexThin == -1)
            continue;
        if(!AcceptHit(silEfficiencies, thin, silIndexThin))
            continue;

        // Energia al llegar a la lamina fina (tras frenar en el gas)
        auto Tat {srim->SlowWithStraggling(srimKey, Tinit, (pointThin - vertex).R())};
        ApplyNaN(Tat);
        if(std::isnan(Tat))
            continue; // se paro en el gas antes de llegar
        double eLossGas {Tinit - Tat};

        auto normal {sils->GetLayer(thin).GetNormal()};
        auto angleWithNormal {TMath::ACos(dirWorldFrame.Unit().Dot(normal.Unit()))};
        silRes->SetParameter(0, sigmaFront);
        auto Tafter {srim->SlowWithStraggling(srimKey + "InSil", Tat, sils->GetLayer(thin).GetUnit().GetThickness(),
                                              angleWithNormal)};
        auto eLossPre {Tat - Tafter};
        auto eLoss {gRandom->Gaus(eLossPre, silRes->Eval(eLossPre))};
        ApplyNaN(eLoss, sils->GetLayer(thin).GetThresholds().at(silIndexThin));
        if(std::isnan(eLoss))
            continue; // no supero el umbral de esta lamina

        res.detected = true;
        res.telescopeCode = t + 1;
        res.layerThin = thin;
        res.silIndexThin = silIndexThin;
        res.point = pointThin;
        res.angle = TMath::ACos(dirWorldFrame.Unit().Dot(normal.Unit())) * TMath::RadToDeg();

        if(Tafter <= 0.)
        {
            // se paro en la lamina fina: energia total = eLoss
            res.punchthrough = false;
            res.energy = eLoss + eLossGas;
        }
        else
        {
            // Punchthrough: busca la lamina gruesa del mismo telescopio
            int silIndexThick {-1};
            ROOT::Math::XYZPoint pointThick;
            std::tie(silIndexThick, pointThick) = sils->FindSPInLayer(thick, vertex, dirWorldFrame);
            res.punchthrough = true;
            res.layerThick = thick;
            if(silIndexThick == -1)
            {
                // atraveso la fina pero no golpea la gruesa: solo se conoce eLoss(fina)
                res.energy = eLoss + eLossGas;
            }
            else
            {
                auto TafterInterGas {srim->SlowWithStraggling(srimKey, Tafter, (pointThin - pointThick).R())};
                auto TafterThick {srim->SlowWithStraggling(srimKey + "InSil", TafterInterGas,
                                                           sils->GetLayer(thick).GetUnit().GetThickness(),
                                                           angleWithNormal)};
                auto eLossBackPre {TafterInterGas - TafterThick};
                auto eLossBack {gRandom->Gaus(eLossBackPre, silRes->Eval(eLossBackPre))};
                ApplyNaN(eLossBack, sils->GetLayer(thick).GetThresholds().at(silIndexThick));
                res.energy = eLoss + eLossGas + (std::isnan(eLossBack) ? 0. : eLossBack);
            }
        }
        return res; // primer telescopio con hit valido
    }
    return res; // detected = false
}

void do_simu_decay(const std::string& beam, const std::string& target, const std::string& light,
                   const std::string& heavy, int neutronPS, int protonPS, double Tbeam, double Ex, bool inspect,
                   bool doAlphaTritonBreakup = false, bool emitNeutron = true, bool useResonance8Li = true,
                   bool useResonance7Li = true, double ExSecondary = 4.63)
{
    // ==================================================================
    // doAlphaTritonBreakup: interruptor maestro. Si es false, el resto de
    //   esta funcion se comporta EXACTAMENTE como antes (otras reacciones:
    //   dd, dt, etc. no se ven afectadas).
    // Si es true, heavy/neutronPS/protonPS pasados por el caller se IGNORAN
    //   para la generacion (se sustituyen por la logica de mas abajo), pero
    //   Ex/heavy se siguen usando tal cual para kinTheo/kin/el umbral, como
    //   ya hacia el codigo original.
    //   emitNeutron     : true  -> canal (d,p): p + 8Li*/7Li* -> ... -> p + n
    //                              (+7Li*) + alfa + t. Usa useResonance8Li
    //                              tal y como antes.
    //                     false -> canal (d,d'): d + 7Li -> d' + 7Li*(Ex) ->
    //                              d' + alfa + t, SIN neutron. En este caso
    //                              useResonance8Li se ignora (no existe nodo
    //                              8Li* posible sin neutron) y solo importa
    //                              useResonance7Li.
    //   useResonance8Li : true  -> se genera primero p + 8Li*(Ex)
    //                     false -> el neutron sale ya en el vertice de
    //                              entrada (nPS=1), heredando Ex como la
    //                              excitacion directa del 7Li*
    //                     (solo aplica si emitNeutron = true)
    //   useResonance7Li : true  -> hay un 7Li* real intermedio antes de
    //                              romper a alfa+t (con excitacion
    //                              ExSecondary si useResonance8Li=true, o
    //                              Ex si useResonance8Li=false o
    //                              emitNeutron=false)
    //                     false -> alfa+t (+n si useResonance8Li=true) se
    //                              generan democraticamente en el mismo paso
    // ==================================================================
    // set batch mode if inspect is false
    if(!inspect)
        gROOT->SetBatch(true);
    // Set whether is PS or not
    bool isPS {(neutronPS > 0) || (protonPS > 0)};
    // Set decayString for files names
    std::string decayStr;
    // Set number of iterations
    const int niter {static_cast<int>(inspect ? 1e6 : (isPS ? 3e7 : 1e7))};
    gRandom->SetSeed(0);
    // Runner: contains utility functions to execute multiple actions as rotate directions
    ActSim::Runner runner(nullptr, nullptr, gRandom, 0);
    // Initialize detectors
    // TPC
    ActRoot::TPCParameters tpc {"Actar"};
    std::cout << "TPC: " << tpc.X() << " " << tpc.Y() << " " << tpc.Z() << '\n';
    // Vertex sampling and beam z variables
    // Each emittance entry: {histogram, fraction of beam time}
    struct EmittanceEntry
    {
        TH3D* hist {};
        double fraction {};
    };
    std::vector<EmittanceEntry> emittances;

    // Give the offsets respect to the front silicon (more stats)
    // Offset defined as Mean silicon point - beam position
    BeamOffsetMap beamOffsets {{"7Li",
                                {
                                    {4.68, 1.0} // fixed beam
                                }}};

    if(beam == "7Li")
    {
        auto f = std::make_unique<TFile>(("../Macros/Emittance/Outputs/histos" + beam + ".root").c_str());
        auto* h = f->Get<TH3D>("h3d");
        if(!h)
            throw std::runtime_error("Could not load beam emittance histogram for 7Li");
        h->SetDirectory(nullptr);
        emittances.push_back({h, 1.0});
    }
    else if(beam == "11Li")
    {
        // Read emittance periods, fractions and beam offsets from .dat file
        std::string datPath {"./Inputs/Efficiencies/beamEmittancePeriods_And_Zoffsets_" + beam + ".dat"};
        std::ifstream finEm(datPath);
        if(!finEm.is_open())
            throw std::runtime_error("Could not open emittance periods file: " + datPath);
        std::string suffix;
        double frac;
        double offset;
        std::vector<BeamOffset>& boVec = beamOffsets[beam]; // create entry for 11Li
        while(finEm >> suffix >> frac >> offset)
        {
            // Load emittance histogram
            std::string fname = "../Macros/Emittance/Outputs/histos" + beam + "_" + suffix + ".root";
            auto f = std::make_unique<TFile>(fname.c_str());
            auto* h = f->Get<TH3D>("h3d");
            if(!h)
                throw std::runtime_error("Could not load beam emittance histogram: " + fname);
            h->SetDirectory(nullptr);
            emittances.push_back({h, frac});
            // Store beam offset for this period
            boVec.push_back({offset, frac});
            std::cout << "Loaded emittance: " << fname << " (fraction = " << frac << ", offset = " << offset << ")\n";
        }
        finEm.close();
    }
    // Silicons
    auto* sils {new ActPhysics::SilSpecs};
    std::string silConfig("silspecs"); // no front silicons, only lateral ones
    sils->ReadFile("../configs/" + silConfig + ".conf");
    // sils->Print();
    const double sigmaSilLat {0.085 / 2.355};   // Si resolution for laterals, around 85 keV FWHM
    const double sigmaSilFront {0.050 / 2.355}; // Si resolution for front, around 50 keV FWHM
    auto silRes = std::make_unique<TF1>(
        "silRes", [=](double* x, double* p) { return p[0] * TMath::Sqrt(x[0] / 5.5); }, 0.0, 100.0, 1);
    std::vector<std::string> silLayers {"l0",
                                        "r0"}; // For the moment only laterals, front not yet available for analysis
    std::vector<std::string> AllsilLayers {"f0", "f1", "f2", "f3", "l0", "r0"};

    std::string filenameSMlat {"../Macros/SilVetos/Outputs/Dists/sms_f0.root"};
    auto fileSMlat {new TFile {filenameSMlat.c_str()}};
    ActPhysics::SilMatrix* smlat =
        fileSMlat->Get<ActPhysics::SilMatrix>("sm5"); // matrix for good distance of left wall
    double silCentreLat = smlat->GetMeanZ({4, 5});
    std::cout << "Silicon left centre at Z = " << silCentreLat << " mm" << std::endl;

    std::string filenameSMfront {"../Macros/SilVetos/Outputs/Dists/sms_f0.root"};
    auto fileSMfront {new TFile {filenameSMfront.c_str()}};
    ActPhysics::SilMatrix* smfront =
        fileSMfront->Get<ActPhysics::SilMatrix>("sm3"); // matrix for good distance of left wall
    double silCentreFront = smfront->GetMeanZ({6, 7, 4});
    std::cout << "Silicon front centre at Z = " << silCentreFront << " mm" << std::endl;

    const double zVertexSigma {0.81}; // From emitance study /always  around the same)

    // We have to centre the silicons with the beam input
    // In real life beam window is not at Z / 2
    // Move lat sils to real placement, I did not do it for f0 because I do not use it for now
    for(auto& [name, layer] : sils->GetLayers())
    {
        if(name == "f0" || name == "f1")
            layer.MoveZTo(silCentreFront, {5});
        if(name == "f2")
            layer.MoveZTo(silCentreFront + 12.5,
                          {0}); // Manually substract 12.5 to get the f2 in the same height as the f3
        if(name == "f3")
            layer.MoveZTo(silCentreFront, {0});
        if(name == "l0" || name == "r0")
            layer.MoveZTo(silCentreLat, {4});
    }
    sils->DrawGeo();
    // Silicon malfunction txt
    std::string silEfficienciesPath {"./Inputs/Efficiencies/silicon_efficiencies_" + beam + ".txt"};
    std::map<std::string, double> silEfficiencies {LoadEfficiencies(silEfficienciesPath)};

    // CUTS ON SILICON ENERGY, depending on particle
    // from the graphical PID cut
    std::string light_name {};
    if(light == "1H")
        light_name = "p";
    else if(light == "2H")
        light_name = "d";
    else if(light == "3H")
        light_name = "t";
    ActRoot::CutsManager<std::string> cuts;
    // read for l0 and r0 and get bigger cut
    std::vector<std::string> silCutLayers {"l0", "r0"};
    std::pair<double, double> eLoss0Cut;
    for(const auto& silLayer : silCutLayers)
    {
        std::string light_name_layer = light_name + silLayer;
        cuts.ReadCut(light_name_layer, TString::Format("../PostAnalysis/Cuts/pid_%s_%s_%s.root", light_name.c_str(),
                                                       silLayer.c_str(), beam.c_str())
                                           .Data());

        if(cuts.GetCut(light_name_layer))
        {
            auto eLoss0Cut_layer = cuts.GetXRange(light_name_layer);
            std::cout << BOLDGREEN << "-> ESil range for " << light_name_layer << " in " << silLayer << ": ["
                      << eLoss0Cut_layer.first << ", " << eLoss0Cut_layer.second << "] MeV" << RESET << '\n';
            if(silLayer == "l0" || eLoss0Cut_layer.second > eLoss0Cut.second)
            {
                eLoss0Cut = eLoss0Cut_layer;
            }
        }
        else
        {
            std::cout << BOLDRED << "Simulation_S2384(): could not read PID cut for " << light << " in " << silLayer
                      << " -> using default eLoss0Cut" << RESET << '\n';
            eLoss0Cut = {0, 1000};
        }
    }

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
    // Tablas SRIM adicionales para el breakup a alfa+triton (nombres fijos,
    // independientes de los strings beam/target/light/heavy pasados)
    std::string alfaName {"4He"};
    std::string tritonName {"3H"};
    if(doAlphaTritonBreakup)
    {
        srim->ReadTable("alfa", path + alfaName + "_" + gas + ".txt");
        srim->ReadTable("alfaInSil", path + alfaName + "_" + silicon + ".txt");
        srim->ReadTable("triton", path + tritonName + "_" + gas + ".txt");
        srim->ReadTable("tritonInSil", path + tritonName + "_" + silicon + ".txt");
    }
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
    auto* kinProton {new ActPhysics::Kinematics {beam, target, "1H", "8Li", Tbeam, Ex}};
    auto* kinDeuteron {new ActPhysics::Kinematics {beam, target, "2H", "7Li", Tbeam, Ex}};

    // ---- Cnstruction of the kinematic generator and, if applicable, the DecayGenerators for the cascade. This is done
    //      ONCE here outside, not in each iteration of the loop. ----
    ActSim::KinematicGenerator* kinGen {nullptr};
    ActSim::DecayGenerator* decayStep2 {nullptr}; // 8Li*->n+7Li* / 8Li*->n+alfa+t / 7Li*->alfa+t
    ActSim::DecayGenerator* decayStep3 {nullptr}; // 7Li*->alfa+t (solo si useResonance8Li && useResonance7Li)
    TGenPhaseSpace* directGenNBody {nullptr};     // solo en los canales democraticos (4 cuerpos con n, 3 sin n)
    int directGenNBodyMult {0};                   // 4 (p,n,alfa,t) o 3 (d',alfa,t)
    TLorentzVector initialLVNBody;

    if(!doAlphaTritonBreakup)
    {
        // Original behaviour
        kinGen = new ActSim::KinematicGenerator {beam, target, light, heavy, protonPS, neutronPS};
    }
    else if(!emitNeutron && useResonance7Li)
    {
        // Canal (d,d'): d + 7Li -> d' + 7Li*(Ex) -> d' + alfa + t (sin neutron)
        kinGen = new ActSim::KinematicGenerator {beam, target, light, "7Li", 0, 0};
        auto li7P {ActPhysics::Particle("7Li")};
        li7P.SetEx(Ex);
        auto d1P {ActPhysics::Particle("4He")};
        auto d2P {ActPhysics::Particle("3H")};
        decayStep2 = new ActSim::DecayGenerator {li7P, d1P, d2P};

        decayStr = Form("dd_7LiEx%.2f", Ex);
    }
    else if(!emitNeutron && !useResonance7Li)
    {
        // Canal (d,d'): d + 7Li -> d' + alfa + t, breakup democratico de 3 cuerpos (sin resonancia real, sin neutron)
        auto* kinAux {new ActPhysics::Kinematics {beam, target, light, "7Li", Tbeam, 0.}};
        auto initLab {kinAux->GetPInitialLab()};
        initialLVNBody = TLorentzVector {initLab.Z(), initLab.Y(), initLab.X(), initLab.E()};
        initialLVNBody *= ActPhysics::Constants::kMeVToGeV;
        double masses3[3] = {ActPhysics::Particle(light).GetMass() * ActPhysics::Constants::kMeVToGeV,
                             ActPhysics::Particle("4He").GetMass() * ActPhysics::Constants::kMeVToGeV,
                             ActPhysics::Particle("3H").GetMass() * ActPhysics::Constants::kMeVToGeV};
        directGenNBody = new TGenPhaseSpace {};
        directGenNBody->SetDecay(initialLVNBody, 3, masses3);
        directGenNBodyMult = 3;

        decayStr = "democratic3Body_dd";
    }
    else if(useResonance8Li && useResonance7Li)
    {
        // p + 8Li*(Ex) -> p + n + 7Li*(ExSecondary) -> p + n + alfa + t
        kinGen = new ActSim::KinematicGenerator {beam, target, light, "8Li", 0, 0};
        auto li8P {ActPhysics::Particle("8Li")};
        li8P.SetEx(Ex);
        auto nP {ActPhysics::Particle("1n")};
        auto li7P {ActPhysics::Particle("7Li")};
        li7P.SetEx(ExSecondary);
        auto d1P {ActPhysics::Particle("4He")};
        auto d2P {ActPhysics::Particle("3H")};
        decayStep2 = new ActSim::DecayGenerator {li8P, nP, li7P};
        decayStep3 = new ActSim::DecayGenerator {li7P, d1P, d2P};

        decayStr = Form("8LiEx%.2f_7LiEx%.2f", Ex, ExSecondary);
    }
    else if(useResonance8Li && !useResonance7Li)
    {
        // p + 8Li*(Ex) -> p + n + alfa + t (breakup democratico del 8Li*)
        kinGen = new ActSim::KinematicGenerator {beam, target, light, "8Li", 0, 0};
        auto li8P {ActPhysics::Particle("8Li")};
        li8P.SetEx(Ex);
        auto nP {ActPhysics::Particle("1n")};
        auto d1P {ActPhysics::Particle("4He")};
        auto d2P {ActPhysics::Particle("3H")};
        decayStep2 = new ActSim::DecayGenerator {li8P, nP, d1P, d2P};

        decayStr = Form("8LiEx%.2f", Ex);
    }
    else if(!useResonance8Li && useResonance7Li)
    {
        // p + n + 7Li*(Ex) -> p + n + alfa + t (n sale ya en el vertice de entrada)
        // OJO: ComputeHeavyMass() resta neutronPS del A del nombre pasado, asi
        // que hay que pasar "8Li" para que, tras restar 1 neutron, el nodo
        // generado sea realmente 7Li. Pasar "7Li" aqui generaria un 6Li.
        kinGen = new ActSim::KinematicGenerator {beam, target, light, "8Li", 0, 1};
        auto li7P {ActPhysics::Particle("7Li")};
        li7P.SetEx(Ex);
        auto d1P {ActPhysics::Particle("4He")};
        auto d2P {ActPhysics::Particle("3H")};
        decayStep2 = new ActSim::DecayGenerator {li7P, d1P, d2P};

        decayStr = Form("7LiEx%.2f", Ex);
    }
    else // !useResonance8Li && !useResonance7Li (con neutron, emitNeutron = true)
    {
        // p + n + alfa + t, breakup democratico a 4 cuerpos, sin resonancias
        // reales. No usa KinematicGenerator: TGenPhaseSpace local directo.
        auto* kinAux {new ActPhysics::Kinematics {beam, target, light, "8Li", Tbeam, 0.}};
        auto initLab {kinAux->GetPInitialLab()};
        initialLVNBody = TLorentzVector {initLab.Z(), initLab.Y(), initLab.X(), initLab.E()};
        initialLVNBody *= ActPhysics::Constants::kMeVToGeV;
        double masses4[4] = {ActPhysics::Particle("1H").GetMass() * ActPhysics::Constants::kMeVToGeV,
                             ActPhysics::Particle("1n").GetMass() * ActPhysics::Constants::kMeVToGeV,
                             ActPhysics::Particle("4He").GetMass() * ActPhysics::Constants::kMeVToGeV,
                             ActPhysics::Particle("3H").GetMass() * ActPhysics::Constants::kMeVToGeV};
        directGenNBody = new TGenPhaseSpace {};
        directGenNBody->SetDecay(initialLVNBody, 4, masses4);
        directGenNBodyMult = 4;

        decayStr = "democratic4Body";
    }

    // cross section sampler
    auto* xs {new ActSim::CrossSection()};
    bool isThereXS {};
    if(neutronPS == 0 && protonPS == 0 && !doAlphaTritonBreakup)
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
    auto hTheta3CMside {HistConfig::ThetaCM.GetHistogram()};
    auto hTheta3CMfront {HistConfig::ThetaCM.GetHistogram()};
    auto hTheta3CMGateHeavy {HistConfig::ThetaCM.GetHistogram()};
    auto hTheta3CMsideGateHeavy {HistConfig::ThetaCM.GetHistogram()};
    auto hThetaCMAll {HistConfig::ThetaCM.GetHistogram()};
    hThetaCMAll->SetTitle("Theta3CM all;#theta_{CM} [#circ];Counts");
    auto hThetaLabAll {HistConfig::ThetaCM.GetHistogram()};
    hThetaLabAll->SetTitle("Theta3Lab all;#theta_{Lab} [#circ];Counts");
    auto hTheta3Lab {HistConfig::ThetaCM.GetHistogram()};
    hTheta3Lab->SetTitle("Theta3Lab;#theta_{Lab} [#circ];Counts");
    auto hTheta3Labside {HistConfig::ThetaCM.GetHistogram()};
    hTheta3Labside->SetTitle("Theta3Lab;#theta_{Lab} [#circ];Counts");
    auto hTheta3Labfront {HistConfig::ThetaCM.GetHistogram()};
    hTheta3Labfront->SetTitle("Theta3Lab;#theta_{Lab} [#circ];Counts");
    auto hTheta3LabGateHeavy {HistConfig::ThetaCM.GetHistogram()};
    hTheta3LabGateHeavy->SetTitle("Theta3Lab;#theta_{Lab} [#circ];Counts");
    auto hTheta3LabsideGateHeavy {HistConfig::ThetaCM.GetHistogram()};
    hTheta3LabsideGateHeavy->SetTitle("Theta3Lab;#theta_{Lab} [#circ];Counts");
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
    auto hExGateHeavy {HistConfig::Ex.GetHistogram()};
    auto hKinRec {HistConfig::Kin.GetHistogram()};
    hKinRec->SetTitle("Reconstructed Kinetic Energy;#theta_{Lab} [#circ];E_{Vertex} [MeV]");
    auto hRP_X {HistConfig::RPx.GetHistogram()};
    hRP_X->SetTitle("Reconstructed RP;X [mm];Counts");
    auto hRP {HistConfig::RP.GetHistogram()};
    hRP->SetTitle("Reconstructed RP;RP [mm];Counts");
    // Debug
    auto hKinDebug {HistConfig::Kin.GetHistogram()};
    hKinDebug->SetTitle("Debug Kinematic Punshthrough;#theta_{Lab} [#circ];E_{Vertex} [MeV]");
    auto* hAlfaEnergyInitial = new TH2D(
        "hAlfaEnergyInitial", "Alfa initial energy;#theta_{Lab} [#circ];E_{Vertex} [MeV]", 180, 0, 180, 150, 0, 50);
    auto* hAlfaEnergyInitialAll = new TH2D(
        "hAlfaEnergyInitialAll", "Alfa initial energy;#theta_{Lab} [#circ];E_{Vertex} [MeV]", 180, 0, 180, 150, 0, 50);
    auto* hAlfaEnergyPost =
        new TH2D("hAlfaEnergyPostSil", "Alfa energy after silicon;#theta_{Lab} [#circ];E_{Vertex} [MeV]", 180, 0, 180,
                 150, 0, 50);
    auto* hEx_asProton = new TH1D("hEx_asProton", "Eex assuming proton;E_{ex} [MeV];Counts", 150, 0, 15);
    auto* hEx_asDeuteron = new TH1D("hEx_asDeuteron", "Eex assuming deuteron;E_{ex} [MeV];Counts", 150, 0, 15);

    // Allow multiple theads
    // std::string tag {""};
    // if(thread > 0)
    //     tag = "_" + std::to_string(thread);

    // File to save data
    TString fileName {TString::Format("./Outputs/%s/Decay/%s_%s_TRIUMF_Eex_%.3f_nPS_%d_pPS_%d_decay_%s.root",
                                      beam.c_str(), target.c_str(), light.c_str(), Ex, neutronPS, protonPS,
                                      decayStr.c_str())};
    // TString fileName {TString::Format("./Outputs/%s/%s_%s_TRIUMF_Eex_%.3f_nPS_%d_pPS_%d%s.root", beam.c_str(),
    //                                   target.c_str(), light.c_str(), Ex, neutronPS, protonPS, tag.c_str())};
    auto outFile {new TFile(fileName, inspect ? "read" : "recreate")};
    auto* outTree {new TTree("SimulationTTree", "A TTree containing only our Eex obtained by simulation")};
    if(inspect)
        outTree->SetDirectory(nullptr);
    double theta3CM_tree {};
    outTree->Branch("theta3CM", &theta3CM_tree);
    double Eex_tree {};
    outTree->Branch("Eex", &Eex_tree);
    double Eex_tree_side {};
    outTree->Branch("Eex_side", &Eex_tree_side);
    double EexGateHeavy_tree {};
    outTree->Branch("EexGateHeavy", &EexGateHeavy_tree);
    double EexGateHeavy_tree_side {};
    outTree->Branch("EexGateHeavy_side", &EexGateHeavy_tree_side);
    double EVertex_tree {};
    outTree->Branch("EVertex", &EVertex_tree);
    double theta3Lab_tree {};
    outTree->Branch("theta3Lab", &theta3Lab_tree);
    double phi3CM_tree {};
    outTree->Branch("phi3CM", &phi3CM_tree);
    double weight_tree {};
    outTree->Branch("weight", &weight_tree);
    double weight_tree_side {};
    outTree->Branch("weight_side", &weight_tree_side);
    double weight_tree_gateHeavy {};
    outTree->Branch("weight_gateHeavy", &weight_tree_gateHeavy);
    double weight_tree_gateHeavy_side {};
    outTree->Branch("weight_gateHeavy_side", &weight_tree_gateHeavy_side);
    // New branches: individual info of detection of alfa and triton (only relevant if doAlphaTritonBreakup=true; remain
    // 0 in other reactions) layerCode: 0 = not detected, 1 = telescope f0/f1, 2 = telescope f2/f3
    int alfaLayerCode_tree {};
    outTree->Branch("alfaLayerCode", &alfaLayerCode_tree);
    double alfaEnergy_tree {};
    outTree->Branch("alfaEnergy", &alfaEnergy_tree);
    double alfaAngle_tree {};
    outTree->Branch("alfaAngle", &alfaAngle_tree);
    int tritonLayerCode_tree {};
    outTree->Branch("tritonLayerCode", &tritonLayerCode_tree);
    double tritonEnergy_tree {};
    outTree->Branch("tritonEnergy", &tritonEnergy_tree);
    double tritonAngle_tree {};
    outTree->Branch("tritonAngle", &tritonAngle_tree);
    // Combined gate: alfa and triton detected (each one in any of its 2 telescopes)
    bool gateAlfaTriton_tree {};
    outTree->Branch("gateAlfaTriton", &gateAlfaTriton_tree);

    // Events stoped inside ACTAR (L1)
    auto* stoppedTree {new TTree("LightMissTree", "Eventos donde el 'light' no llega a reconstruirse")};
    if(inspect)
        stoppedTree->SetDirectory(nullptr);
    int lightStatus_treeStop {}; // 0 = ningun silicio en la trayectoria, 1 = se para en el gas antes de llegar
    stoppedTree->Branch("lightStatus", &lightStatus_treeStop);
    double lightTheta3Lab_treeStop {};
    stoppedTree->Branch("theta3Lab", &lightTheta3Lab_treeStop);
    double lightT3Lab_treeStop {};
    stoppedTree->Branch("T3Lab", &lightT3Lab_treeStop);
    double lightRangeInGas_treeStop {};
    stoppedTree->Branch("rangeInGas", &lightRangeInGas_treeStop);
    double lightEx_treeStop {};
    stoppedTree->Branch("Eex", &lightEx_treeStop);
    // reutiliza las variables de alfa/triton (ya calculadas arriba, antes de estos continues)
    stoppedTree->Branch("alfaLayerCode", &alfaLayerCode_tree);
    stoppedTree->Branch("alfaEnergy", &alfaEnergy_tree);
    stoppedTree->Branch("alfaAngle", &alfaAngle_tree);
    double alfaEnergyTruth_treeStop {};
    stoppedTree->Branch("alfaEnergyTruth", &alfaEnergyTruth_treeStop);
    double alfaThetaLabTruth_treeStop {};
    stoppedTree->Branch("alfaAngleTruth", &alfaThetaLabTruth_treeStop);
    stoppedTree->Branch("tritonLayerCode", &tritonLayerCode_tree);
    stoppedTree->Branch("tritonEnergy", &tritonEnergy_tree);
    stoppedTree->Branch("tritonAngle", &tritonAngle_tree);
    double tritonEnergyTruth_treeStop {};
    stoppedTree->Branch("tritonEnergyTruth", &tritonEnergyTruth_treeStop);
    double tritonThetaLabTruth_treeStop {};
    stoppedTree->Branch("tritonAngleTruth", &tritonThetaLabTruth_treeStop);
    stoppedTree->Branch("gateAlfaTriton", &gateAlfaTriton_tree);

    // Counter of events in which the proton did not have enough energy to exit the ACTAR gas volume


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
        double beamOffset = PickBeamOffset(beam, beamOffsets);
        double zVertexMeanEvt = silCentreFront - beamOffset;
        // Pick emittance histogram according to beam-time fractions
        TH3D* hBeam = [&]() -> TH3D*
        {
            double rEm = gRandom->Uniform();
            double accEm = 0.0;
            for(auto& em : emittances)
            {
                accEm += em.fraction;
                if(rEm < accEm)
                    return em.hist;
            }
            return emittances.back().hist;
        }();
        auto [start, vertex] {SampleVertex(zVertexMeanEvt, zVertexSigma, hBeam, tpc.X())};
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
        // Productos del breakup alfa+triton (y neutron, invisible en Si).
        // Se rellenan mas abajo solo si doAlphaTritonBreakup=true.
        TLorentzVector alfaLV {}, tritonLV {}, neutronLV {};
        bool haveAlphaTriton {false};
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
        else if(directGenNBody)
        {
            // Breakup democratico (p+n+alfa+t si emitNeutron, o d'+alfa+t si
            // !emitNeutron), sin resonancias reales: T3Lab/theta3Lab/phi3Lab
            // TAMBIEN salen de aqui, no hay kinGen en absoluto en este modo.
            weight = directGenNBody->Generate();
            TLorentzVector lightLV {*directGenNBody->GetDecay(0)}; // p (con n) o d' (sin n)
            int idxAlfa {1};
            int idxTriton {2};
            if(directGenNBodyMult == 4)
            {
                neutronLV = *directGenNBody->GetDecay(1);
                idxAlfa = 2;
                idxTriton = 3;
            }
            alfaLV = *directGenNBody->GetDecay(idxAlfa);
            tritonLV = *directGenNBody->GetDecay(idxTriton);
            // TGenPhaseSpace trabaja en GeV con las masas que le dimos; reescalar a MeV
            lightLV *= 1. / ActPhysics::Constants::kMeVToGeV;
            if(directGenNBodyMult == 4)
                neutronLV *= 1. / ActPhysics::Constants::kMeVToGeV;
            alfaLV *= 1. / ActPhysics::Constants::kMeVToGeV;
            tritonLV *= 1. / ActPhysics::Constants::kMeVToGeV;
            haveAlphaTriton = true;

            theta3Lab = lightLV.Theta();
            phi3Lab = lightLV.Phi();
            T3Lab = lightLV.E() - lightLV.M();
            theta3LabSampled = theta3Lab;
            ApplyThetaRes(theta3Lab);

            // OJO: 'kin' aqui es el binario original (heavy/Ex pasados a la
            // funcion), que no describe exactamente este breakup democratico
            // (no hay 2-cuerpos real). Se usa igualmente como referencia
            // aproximada para theta3CM, igual que el resto del pipeline.
            theta3CMBefore = kin->ReconstructTheta3CMFromLab(T3Lab, theta3LabSampled) * TMath::RadToDeg();
            theta3CM = kin->ReconstructTheta3CMFromLab(T3Lab, theta3Lab);
            phi3CM = phi3Lab;

            // Ya no existe un unico "heavy": se deja a 0, no se usa mas abajo
            theta4Lab = 0.;
            phi4Lab = 0.;
            T4Lab = 0.;
        }
        else
        {
            // Sample kinematics generator
            kinGen->SetBeamAndExEnergies(TbeamCorr, randEx);
            weight = kinGen->Generate();
            if(neutronPS == 0 && protonPS == 0 && !doAlphaTritonBreakup)
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

            // Heavy (nodo intermedio: 8Li*, 7Li*, o "heavy" original si no hay breakup)
            theta4Lab = LorenztVector4->Theta();
            phi4Lab = LorenztVector4->Phi();
            T4Lab = LorenztVector4->E() - LorenztVector4->M();

            // ---- Cascada de decaimiento del nodo "heavy" a productos detectables ----
            if(doAlphaTritonBreakup)
            {
                if(!emitNeutron && useResonance7Li)
                {
                    // Canal (d,d'): 7Li* -> alfa + t directamente, sin neutron
                    decayStep2->SetDecay(T4Lab, theta4Lab, phi4Lab);
                    double w2 {decayStep2->Generate()};
                    alfaLV = *decayStep2->GetLorentzVector(0);
                    tritonLV = *decayStep2->GetLorentzVector(1);
                    weight *= w2;
                }
                else if(useResonance8Li && useResonance7Li)
                {
                    // 8Li* -> n + 7Li*
                    decayStep2->SetDecay(T4Lab, theta4Lab, phi4Lab);
                    double w2 {decayStep2->Generate()};
                    neutronLV = *decayStep2->GetLorentzVector(0);
                    auto* li7starLV {decayStep2->GetLorentzVector(1)};
                    // 7Li* -> alfa + t
                    double T7 {li7starLV->E() - li7starLV->M()};
                    decayStep3->SetDecay(T7, li7starLV->Theta(), li7starLV->Phi());
                    double w3 {decayStep3->Generate()};
                    alfaLV = *decayStep3->GetLorentzVector(0);
                    tritonLV = *decayStep3->GetLorentzVector(1);
                    weight *= w2 * w3;
                }
                else if(useResonance8Li && !useResonance7Li)
                {
                    // 8Li* -> n + alfa + t (breakup democratico de 3 cuerpos)
                    decayStep2->SetDecay(T4Lab, theta4Lab, phi4Lab);
                    double w2 {decayStep2->Generate()};
                    neutronLV = *decayStep2->GetLorentzVector(0);
                    alfaLV = *decayStep2->GetLorentzVector(1);
                    tritonLV = *decayStep2->GetLorentzVector(2);
                    weight *= w2;
                }
                else // !useResonance8Li && useResonance7Li (emitNeutron = true)
                {
                    // el neutron ya salio directo de kinGen (3 cuerpos: p,n,7Li*)
                    neutronLV = *kinGen->GetLorentzVector(2);
                    // 7Li* -> alfa + t
                    decayStep2->SetDecay(T4Lab, theta4Lab, phi4Lab);
                    double w2 {decayStep2->Generate()};
                    alfaLV = *decayStep2->GetLorentzVector(0);
                    tritonLV = *decayStep2->GetLorentzVector(1);
                    weight *= w2;
                }
                haveAlphaTriton = true;
            }
        }
        double alfaThetaLabTruth {alfaLV.Theta() * TMath::RadToDeg()};
        double alfaEnergyTruth {alfaLV.E() - alfaLV.M()};
        double tritonThetaLabTruth {tritonLV.Theta() * TMath::RadToDeg()};
        double tritonEnergyTruth {tritonLV.E() - tritonLV.M()};
        // Fill debug alfa energy
        hAlfaEnergyInitialAll->Fill(alfaLV.Theta() * TMath::RadToDeg(), alfaLV.E() - alfaLV.M());
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

        // Threshold L1, particles that stop in actar. Check before doing the continues
        double rangeInGas {srim->EvalRange("light", T3Lab)};
        ROOT::Math::XYZPoint finalPointGas {vertex + rangeInGas * dirWorldFrame.Unit()};
        if(0 <= finalPointGas.X() && finalPointGas.X() <= 256 && 0 <= finalPointGas.Y() && finalPointGas.Y() <= 256 &&
           std::abs(vertex.Z() - finalPointGas.Z()) < (235. / 2.))
        {

        } // the z axis is referred to experimental values, it is not used for L1

        // ---- Tracking del alfa y del triton (adelantado antes de los 'continue' del light) ----
        // Se calcula aqui, independientemente de lo que le pase al 'light', para
        // no perder esta informacion en los eventos en los que el 'light' no
        // llega a ningun silicio o se para en el gas.
        alfaLayerCode_tree = 0;
        alfaEnergy_tree = -1000;
        alfaAngle_tree = -1000;
        tritonLayerCode_tree = 0;
        tritonEnergy_tree = -1000;
        tritonAngle_tree = -1000;
        gateAlfaTriton_tree = false;
        bool silIndexHeavy_equivalent {false}; // true si pasa el gate combinado
        // Info "verdad" (Lorentz vector), independiente de si detectan en el Si o no
        double alfaEnergyTruth_tree {-1000};
        double alfaAngleTruth_tree {-1000};
        double tritonEnergyTruth_tree {-1000};
        double tritonAngleTruth_tree {-1000};
        if(doAlphaTritonBreakup && haveAlphaTriton)
        {
            alfaEnergyTruth_tree = alfaLV.E() - alfaLV.M();
            alfaAngleTruth_tree = alfaLV.Theta() * TMath::RadToDeg();
            tritonEnergyTruth_tree = tritonLV.E() - tritonLV.M();
            tritonAngleTruth_tree = tritonLV.Theta() * TMath::RadToDeg();

            ROOT::Math::XYZVector alfaBeamFrame {TMath::Cos(alfaLV.Theta()),
                                                 TMath::Sin(alfaLV.Theta()) * TMath::Sin(alfaLV.Phi()),
                                                 TMath::Sin(alfaLV.Theta()) * TMath::Cos(alfaLV.Phi())};
            ROOT::Math::XYZVector tritonBeamFrame {TMath::Cos(tritonLV.Theta()),
                                                   TMath::Sin(tritonLV.Theta()) * TMath::Sin(tritonLV.Phi()),
                                                   TMath::Sin(tritonLV.Theta()) * TMath::Cos(tritonLV.Phi())};
            auto alfaWorldFrame {runner.RotateToWorldFrame(alfaBeamFrame, beamDir)};
            auto tritonWorldFrame {runner.RotateToWorldFrame(tritonBeamFrame, beamDir)};

            double TalfaLab {alfaLV.E() - alfaLV.M()};
            double TtritonLab {tritonLV.E() - tritonLV.M()};

            auto alfaHit {TrackFrontWall("alfa", TalfaLab, alfaWorldFrame, vertex, sils, srim, silRes.get(),
                                         sigmaSilFront, silEfficiencies)};
            auto tritonHit {TrackFrontWall("triton", TtritonLab, tritonWorldFrame, vertex, sils, srim, silRes.get(),
                                           sigmaSilFront, silEfficiencies)};

            alfaLayerCode_tree = alfaHit.telescopeCode;
            alfaEnergy_tree = alfaHit.detected ? alfaHit.energy : -1000;
            alfaAngle_tree = alfaHit.detected ? alfaHit.angle : -1000;
            tritonLayerCode_tree = tritonHit.telescopeCode;
            tritonEnergy_tree = tritonHit.detected ? tritonHit.energy : -1000;
            tritonAngle_tree = tritonHit.detected ? tritonHit.angle : -1000;
            gateAlfaTriton_tree = alfaHit.detected && tritonHit.detected;
            silIndexHeavy_equivalent = gateAlfaTriton_tree;
        }

        // How to check whether tracks would read the silicons with new class:
        int silIndex0 = -1;
        ROOT::Math::XYZPoint silPoint0;
        std::string layer0;
        for(auto layer : silLayers)
        {
            std::tie(silIndex0, silPoint0) = sils->FindSPInLayer(layer, vertex, dirWorldFrame);
            if(silIndex0 != -1)
            {
                layer0 = layer;
                break;
            }
        }
        // Fix silion resolution depending on hte layer hit
        silRes->SetParameter(0, layer0 == "f0" ? sigmaSilFront : sigmaSilLat);
        if(silIndex0 == -1)
        {
            if(0 <= finalPointGas.X() && finalPointGas.X() <= 256 && 0 <= finalPointGas.Y() &&
               finalPointGas.Y() <= 256 && std::abs(vertex.Z() - finalPointGas.Z()) < (235. / 2.))
            {
                // El 'light' no apunta a ningun silicio: se guarda en el arbol de
                // eventos "perdidos" con la info de alfa/triton ya calculada arriba.
                lightStatus_treeStop = 0;
                lightTheta3Lab_treeStop = theta3Lab * TMath::RadToDeg();
                lightT3Lab_treeStop = T3Lab;
                lightRangeInGas_treeStop = rangeInGas;
                lightEx_treeStop = kin->ReconstructExcitationEnergy(T3Lab, theta3Lab);
                alfaEnergyTruth_treeStop = alfaEnergyTruth_tree;
                alfaThetaLabTruth_treeStop = alfaAngleTruth_tree;
                tritonEnergyTruth_treeStop = tritonEnergyTruth_tree;
                tritonThetaLabTruth_treeStop = tritonAngleTruth_tree;
                stoppedTree->Fill();
            }
            continue;
        }
        // Slow down light in gas
        auto T3AtSil {srim->SlowWithStraggling("light", T3Lab, (silPoint0 - vertex).R())};
        // Check if stopped
        ApplyNaN(T3AtSil);
        if(std::isnan(T3AtSil))
        {
            if(0 <= finalPointGas.X() && finalPointGas.X() <= 256 && 0 <= finalPointGas.Y() &&
               finalPointGas.Y() <= 256 && std::abs(vertex.Z() - finalPointGas.Z()) < (235. / 2.))
            {
                // El 'light' se para en el gas: se guarda en el arbol de eventos
                // "perdidos" con la info de alfa/triton ya calculada arriba.
                lightStatus_treeStop = 1;
                lightTheta3Lab_treeStop = theta3Lab * TMath::RadToDeg();
                lightT3Lab_treeStop = T3Lab;
                lightRangeInGas_treeStop = rangeInGas;
                lightEx_treeStop = kin->ReconstructExcitationEnergy(T3Lab, theta3Lab);
                alfaEnergyTruth_treeStop = alfaEnergyTruth_tree;
                alfaThetaLabTruth_treeStop = alfaAngleTruth_tree;
                tritonEnergyTruth_treeStop = tritonEnergyTruth_tree;
                tritonThetaLabTruth_treeStop = tritonAngleTruth_tree;
                stoppedTree->Fill();
            }
            continue;
        }
        // Slow down in silicon
        auto normal {sils->GetLayer(layer0).GetNormal()};
        auto angleWithNormal {TMath::ACos(dirWorldFrame.Unit().Dot(normal.Unit()))};
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
            std::tie(silIndex1, silPoint1) = sils->FindSPInLayer("f1", vertex, dirWorldFrame);
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
        if(silIndexHeavy_equivalent)
        {
            hTheta3CMGateHeavy->Fill(theta3CMBefore);
            hTheta3LabGateHeavy->Fill(theta3Lab * TMath::RadToDeg());
            if(layer0 == "l0" || layer0 == "r0")
            {
                hTheta3CMsideGateHeavy->Fill(theta3CMBefore);
                hTheta3LabsideGateHeavy->Fill(theta3Lab * TMath::RadToDeg());
            }
        }
        else
        {
            EexGateHeavy_tree = -1000; // if gate not passed, fill with -1000 to be able to do gate in analysis
            EexGateHeavy_tree_side = -1000;
            weight_tree_gateHeavy = -1000;
            weight_tree_gateHeavy_side = -1000;
        }
        // Check if corresponds to a hit when the detector was off or on
        if(!AcceptHit(silEfficiencies, layer0, silIndex0))
        {
            continue; // if not accepted, go to next iteration
        }
        // Reconstruct!
        if(T3AfterSil0 > 0)
        {
            hKinDebug->Fill(theta3Lab * TMath::RadToDeg(), T3Lab);
        }
        bool isOnlyFirstWall {(T3AfterSil0 == 0)}; // Only analyse the first wall, no punchthrough, as in the experiment
        bool isOk {(T3AfterSil0 == 0 || T3AfterSil1 == 0)}; // no punchthrouhg
        bool cutELoss0 {true};                              // for f0 not yet implemented the graphical cuts
        if(layer0 == "r0" || layer0 == "l0")
            cutELoss0 = (eLoss0Cut.first <= eLoss0 && eLoss0 <= eLoss0Cut.second); // graphical cuts on experimental PID
        if(isOnlyFirstWall && cutELoss0)
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
            hEx->Fill(ExRec, weight);                            // To get real counts weigth * alpha
            hRP_X->Fill(vertex.X());
            hRP->Fill(vertex.X(), vertex.Y());
            // Debug proton/deuteron contamination
            hEx_asDeuteron->Fill(kinDeuteron->ReconstructExcitationEnergy(T3Rec, theta3Lab), weight);
            hEx_asProton->Fill(kinProton->ReconstructExcitationEnergy(T3Rec, theta3Lab), weight);
            if(layer0 == "f0")
            {
                hSPf0->Fill(silPoint0.Y(), silPoint0.Z());
                hTheta3CMfront->Fill(theta3CMBefore);
                hTheta3Labfront->Fill(theta3Lab * TMath::RadToDeg());
            }
            if(layer0 == "l0")
            {
                hSPl0->Fill(silPoint0.X(), silPoint0.Z());
                hTheta3CMside->Fill(theta3CMBefore);
                hTheta3Labside->Fill(theta3Lab * TMath::RadToDeg());
            }
            if(layer0 == "r0")
            {
                hSPr0->Fill(silPoint0.X(), silPoint0.Z());
                hTheta3CMside->Fill(theta3CMBefore);
                hTheta3Labside->Fill(theta3Lab * TMath::RadToDeg());
            }
            hTheta3CM->Fill(theta3CMBefore); // only thetaCm that enter our cuts
            hTheta3Lab->Fill(theta3Lab * TMath::RadToDeg());
            // write to TTree
            Eex_tree = ExRec;
            theta3CM_tree = theta3CM * TMath::RadToDeg();
            EVertex_tree = T3Rec;
            theta3Lab_tree = theta3Lab * TMath::RadToDeg();
            phi3CM_tree = phi3CM;
            weight_tree = weight;
            if(layer0 == "l0" || layer0 == "r0")
            {
                Eex_tree_side = ExRec;
                weight_tree_side = weight;
                if(silIndexHeavy_equivalent)
                {
                    EexGateHeavy_tree_side = ExRec;
                    weight_tree_gateHeavy_side = weight;
                }
            }
            else
            {
                Eex_tree_side = -1000; // if not side, fill with -1000 to be able to do gate in analysis
                weight_tree_side = -1000;
                EexGateHeavy_tree_side = -1000;
                weight_tree_gateHeavy_side = -1000;
            }
            if(silIndexHeavy_equivalent)
            {
                EexGateHeavy_tree = ExRec;
                weight_tree_gateHeavy = weight;
            }
            outTree->Fill();
        }
        hAlfaEnergyInitial->Fill(alfaAngle_tree, alfaLV.E() - alfaLV.M());
        hAlfaEnergyPost->Fill(alfaAngle_tree, alfaEnergy_tree);
    }


    // Compute efficiency side, front, heavy gate and total
    auto* effCM {new TEfficiency {*hTheta3CM, *hThetaCMAll}};
    effCM->SetNameTitle("effCM", " #epsilon_{TOT} (#theta_{CM});#epsilon;#theta_{CM} [#circ]");
    auto* effLab {new TEfficiency {*hTheta3Lab, *hThetaLabAll}};
    effLab->SetNameTitle("effLab", "#epsilon_{TOT} (#theta_{Lab});#epsilon;#theta_{Lab} [#circ]");

    auto* effCMside {new TEfficiency {*hTheta3CMside, *hThetaCMAll}};
    effCMside->SetNameTitle("effCMside", " #epsilon_{side} (#theta_{CM});#epsilon;#theta_{CM} [#circ]");
    auto* effLabside {new TEfficiency {*hTheta3Labside, *hThetaLabAll}};
    effLabside->SetNameTitle("effLabside", "#epsilon_{side} (#theta_{Lab});#epsilon;#theta_{Lab} [#circ]");

    auto* effCMfront {new TEfficiency {*hTheta3CMfront, *hThetaCMAll}};
    effCMfront->SetNameTitle("effCMfront", " #epsilon_{front} (#theta_{CM});#epsilon;#theta_{CM} [#circ]");
    auto* effLabfront {new TEfficiency {*hTheta3Labfront, *hThetaLabAll}};
    effLabfront->SetNameTitle("effLabfront", "#epsilon_{front} (#theta_{Lab});#epsilon;#theta_{Lab} [#circ]");

    auto* effCMgateHeavy {new TEfficiency {*hTheta3CMGateHeavy, *hThetaCMAll}};
    effCMgateHeavy->SetNameTitle("effCMgateHeavy", " #epsilon_{heavy gate} (#theta_{CM});#epsilon;#theta_{CM} [#circ]");
    auto* effLabgateHeavy {new TEfficiency {*hTheta3LabGateHeavy, *hThetaLabAll}};
    effLabgateHeavy->SetNameTitle("effLabgateHeavy",
                                  "#epsilon_{heavy gate} (#theta_{Lab});#epsilon;#theta_{Lab} [#circ]");
    auto* effCMsideGateHeavy {new TEfficiency {*hTheta3CMsideGateHeavy, *hThetaCMAll}};
    effCMsideGateHeavy->SetNameTitle("effCMsideGateHeavy",
                                     " #epsilon_{side} (#theta_{CM});#epsilon;#theta_{CM} [#circ]");
    auto* effLabsideGateHeavy {new TEfficiency {*hTheta3LabsideGateHeavy, *hThetaLabAll}};
    effLabsideGateHeavy->SetNameTitle("effLabsideGateHeavy",
                                      "#epsilon_{side} (#theta_{Lab});#epsilon;#theta_{Lab} [#circ]");

    // Reportar (siempre) el numero de protones que se quedaron sin energia
    // suficiente para salir del volumen de gas de ACTAR

    // SAVING
    if(!inspect)
    {
        outFile->cd();
        outTree->Write();
        stoppedTree->Write();
        effCM->Write();
        effLab->Write();
        effCMside->Write();
        effLabside->Write();
        effCMfront->Write();
        effLabfront->Write();
        effCMgateHeavy->Write();
        effLabgateHeavy->Write();
        hRP->Write("hRP");
        outFile->Close();
        delete outFile;
        outFile = nullptr;
    }
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
        c1->cd(3);
        hKinDebug->DrawClone("colz");

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

        auto* cDebugAlfa {new TCanvas {"cDebugAlfa", "Debug Alfa"}};
        cDebugAlfa->DivideSquare(3);
        cDebugAlfa->cd(1);
        hAlfaEnergyInitialAll->DrawClone("colz");
        cDebugAlfa->cd(2);
        hAlfaEnergyInitial->DrawClone("colz");
        cDebugAlfa->cd(3);
        hAlfaEnergyPost->DrawClone("colz");

        auto* cDebugEx {new TCanvas {"cDebugEx", "Debug Ex"}};
        cDebugEx->DivideSquare(2);
        cDebugEx->cd(1);
        hEx_asDeuteron->DrawClone("hist");
        cDebugEx->cd(2);
        hEx_asProton->DrawClone("hist");
    }

    // deleting news
    delete srim;
    if(isThereXS)
        delete xs;
}
#endif