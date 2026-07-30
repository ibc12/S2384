#include "ActColors.h"
#include "ActConstants.h"
#include "ActDecayGenerator.h"
#include "ActKinematicGenerator.h"
#include "ActKinematics.h"
#include "ActParticle.h"

#include "TCanvas.h"
#include "TGenPhaseSpace.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TMath.h"

#include <iostream>
#include <map>
#include <string>

// ======================================================================
// Los 4 escenarios fisicos posibles se derivan de 2 preguntas binarias
// independientes; no hace falta un enum, con 2 booleanos basta:
//
//  kUseResonance8Li = true, kUseResonance7Li = true
//      7Li+d -> p + 8Li*(Ex8) -> p + n + 7Li*(Ex7) -> p + n + alfa + t
//      (dos resonancias reales, cascada de 3 pasos)
//
//  kUseResonance8Li = true, kUseResonance7Li = false
//      7Li+d -> p + 8Li*(Ex8) -> p + n + alfa + t
//      (8Li real, pero rompe democraticamente a n+alfa+t sin pasar por un 7Li*)
//
//  kUseResonance8Li = false, kUseResonance7Li = true
//      7Li+d -> p + n + 7Li*(Ex7) -> p + n + alfa + t
//      (el neutron sale ya en el vertice de entrada -> nPS=1 en kinGen;
//       el 7Li* SI es una resonancia real que luego rompe a alfa+t)
//
//  kUseResonance8Li = false, kUseResonance7Li = false
//      7Li+d -> p + n + alfa + t
//      (breakup democratico a 4 cuerpos, sin ninguna resonancia intermedia;
//       no pasa por KinematicGenerator en absoluto, porque esa clase solo
//       sabe anadir n/p sueltos, no alfa/t)
// ======================================================================

// Histogramas theta/energia para una particula generica
struct ParticleHistos
{
    TH1D* theta;
    TH1D* energy;
    TH2D* kin;
    ParticleHistos(const std::string& name, double emax)
    {
        theta = new TH1D(("hTheta_" + name).c_str(), (name + " Lab Angle;#theta_{Lab} [#circ];Counts").c_str(), 100, 0,
                         180);
        energy =
            new TH1D(("hEnergy_" + name).c_str(), (name + " Lab Energy;E_{Lab} [MeV];Counts").c_str(), 100, 0, emax);
        kin = new TH2D(("hKin_" + name).c_str(), (name + " Lab Kinematics;#theta_{Lab} [#circ];E_{Lab} [MeV]").c_str(),
                       100, 0, 180, 100, 0, emax);
    }
    void Fill(const TLorentzVector& lv, double weight)
    {
        theta->Fill(lv.Theta() * TMath::RadToDeg(), weight);
        energy->Fill(lv.E() - lv.M(), weight);
        kin->Fill(lv.Theta() * TMath::RadToDeg(), lv.E() - lv.M(), weight);
    }
};

// Construye a mano el 4-vector inicial beam+target, con el mismo criterio de
// ejes y unidades que usa ActKinematicGenerator::GetInitialState() por dentro
// (necesario solo para kDirect_4Body, que no pasa por KinematicGenerator)
TLorentzVector BuildInitialLV(ActPhysics::Kinematics* kin)
{
    auto init {kin->GetPInitialLab()};
    TLorentzVector lv {init.Z(), init.Y(), init.X(), init.E()};
    lv *= ActPhysics::Constants::kMeVToGeV;
    return lv;
}

void testDecayAlfaTriton_generalized()
{
    // ==================================================================
    // ---- CONFIGURACION: edita aqui, no hace falta pasar argumentos ----
    // ==================================================================
    bool kUseResonance8Li {false}; // true: 8Li real (p+8Li*, 2 cuerpos) primero
                                   // false: el neutron sale ya en el vertice de entrada (nPS=1)
    bool kUseResonance7Li {true}; // true: 7Li* real intermedio antes de romper a alfa+t
                                   // false: alfa+t se generan junto con n en el mismo breakup

    double Ex8Li {7.1};  // Energia de excitacion del 8Li (solo se usa si kUseResonance8Li=true)
    double Ex7Li {4.63}; // Energia de excitacion del 7Li (solo se usa si kUseResonance7Li=true)
    // ==================================================================

    // Particulas y energia de haz, comunes a todos los modos
    std::string beam {"7Li"};
    std::string target {"2H"};
    std::string light {"1H"};
    double TBeam {7 * 7.558}; // MeV
    const int niter {100000};

    // Histogramas genericos por particula
    std::map<std::string, ParticleHistos*> h;
    h["p"] = new ParticleHistos("p", 20);
    h["n"] = new ParticleHistos("n", 20);
    h["alfa"] = new ParticleHistos("alfa", 50);
    h["t"] = new ParticleHistos("t", 50);

    // Chequeo universal de conservacion de energia total: debe dar 0 en LOS 4 MODOS.
    // A diferencia de "EnergyDifference" (que solo mide el Q-valor de un paso concreto),
    // este balance valida la cadena COMPLETA de generacion.
    auto* hBalance = new TH1D("hBalance", "Balance E_{inicial} - #sum E_{final};#Delta E [MeV];Counts", 200, -1, 1);

    double Einitial {TBeam + ActPhysics::Particle(beam).GetMass() + ActPhysics::Particle(target).GetMass()};

    // ---- Construccion de generadores segun el modo (UNA VEZ, fuera del loop) ----
    ActSim::KinematicGenerator* kinGen {nullptr};
    ActSim::DecayGenerator* decayStep2 {nullptr}; // 8Li*->n+7Li*  o  8Li*->n+alfa+t  o  7Li*->alfa+t
    ActSim::DecayGenerator* decayStep3 {nullptr}; // 7Li*->alfa+t (solo en kCascade_8Li_7Li)
    TGenPhaseSpace* directGen4Body {nullptr};
    TLorentzVector initialLV4Body;

    if(kUseResonance8Li && kUseResonance7Li)
    {
        // p + 8Li*(Ex8) -> p + n + 7Li*(Ex7) -> p + n + alfa + t
        kinGen = new ActSim::KinematicGenerator {beam, target, light, "8Li", 0, 0};
        auto li8P {ActPhysics::Particle("8Li")};
        li8P.SetEx(Ex8Li);
        auto nP {ActPhysics::Particle("1n")};
        auto li7P {ActPhysics::Particle("7Li")};
        li7P.SetEx(Ex7Li);
        auto d1P {ActPhysics::Particle("4He")};
        auto d2P {ActPhysics::Particle("3H")};
        decayStep2 = new ActSim::DecayGenerator {li8P, nP, li7P}; // 8Li* -> n + 7Li*
        decayStep3 = new ActSim::DecayGenerator {li7P, d1P, d2P}; // 7Li* -> alfa + t
    }
    else if(kUseResonance8Li && !kUseResonance7Li)
    {
        // p + 8Li*(Ex8) -> p + n + alfa + t (breakup democratico del 8Li*)
        kinGen = new ActSim::KinematicGenerator {beam, target, light, "8Li", 0, 0};
        auto li8P {ActPhysics::Particle("8Li")};
        li8P.SetEx(Ex8Li);
        auto nP {ActPhysics::Particle("1n")};
        auto d1P {ActPhysics::Particle("4He")};
        auto d2P {ActPhysics::Particle("3H")};
        decayStep2 = new ActSim::DecayGenerator {li8P, nP, d1P, d2P}; // 8Li* -> n+alfa+t
    }
    else if(!kUseResonance8Li && kUseResonance7Li)
    {
        // p + n + 7Li*(Ex7) -> p + n + alfa + t (n sale ya en el vertice de entrada)
        // OJO: ComputeHeavyMass() resta nPS del A del nombre pasado, asi que hay
        // que pasar "8Li" (A=8) para que, tras restar 1 neutron, el nodo generado
        // sea realmente 7Li (A=7). Pasar "7Li" aqui generaria por error un 6Li.
        kinGen = new ActSim::KinematicGenerator {beam, target, light, "8Li", 0, 1}; // nPS=1
        auto li7P {ActPhysics::Particle("7Li")};
        li7P.SetEx(Ex7Li);
        auto d1P {ActPhysics::Particle("4He")};
        auto d2P {ActPhysics::Particle("3H")};
        decayStep2 = new ActSim::DecayGenerator {li7P, d1P, d2P}; // 7Li* -> alfa + t
    }
    else // !kUseResonance8Li && !kUseResonance7Li
    {
        // p + n + alfa + t, breakup democratico a 4 cuerpos, sin resonancias.
        // No usa KinematicGenerator: se genera todo con un TGenPhaseSpace local
        auto* kinAux {new ActPhysics::Kinematics {beam, target, light, "8Li", TBeam, 0.}};
        initialLV4Body = BuildInitialLV(kinAux);
        double masses[4] = {ActPhysics::Particle("1H").GetMass() * ActPhysics::Constants::kMeVToGeV,
                            ActPhysics::Particle("1n").GetMass() * ActPhysics::Constants::kMeVToGeV,
                            ActPhysics::Particle("4He").GetMass() * ActPhysics::Constants::kMeVToGeV,
                            ActPhysics::Particle("3H").GetMass() * ActPhysics::Constants::kMeVToGeV};
        directGen4Body = new TGenPhaseSpace {};
        directGen4Body->SetDecay(initialLV4Body, 4, masses);
    }

    // ---- Loop de generacion ----
    for(int i = 0; i < niter; ++i)
    {
        double weight {1.};
        TLorentzVector pLV, nLV, alfaLV, tritonLV;

        if(kUseResonance8Li && kUseResonance7Li)
        {
            kinGen->SetBeamAndExEnergies(TBeam, Ex8Li);
            double w1 {kinGen->Generate()};
            pLV = *kinGen->GetLorentzVector(0);
            auto* li8starLV {kinGen->GetLorentzVector(1)};

            double T4 {li8starLV->E() - li8starLV->M()};
            decayStep2->SetDecay(T4, li8starLV->Theta(), li8starLV->Phi());
            double w2 {decayStep2->Generate()};
            nLV = *decayStep2->GetLorentzVector(0);
            auto* li7starLV {decayStep2->GetLorentzVector(1)};

            double T7 {li7starLV->E() - li7starLV->M()};
            decayStep3->SetDecay(T7, li7starLV->Theta(), li7starLV->Phi());
            double w3 {decayStep3->Generate()};
            alfaLV = *decayStep3->GetLorentzVector(0);
            tritonLV = *decayStep3->GetLorentzVector(1);

            weight = w1 * w2 * w3;
        }
        else if(kUseResonance8Li && !kUseResonance7Li)
        {
            kinGen->SetBeamAndExEnergies(TBeam, Ex8Li);
            double w1 {kinGen->Generate()};
            pLV = *kinGen->GetLorentzVector(0);
            auto* li8starLV {kinGen->GetLorentzVector(1)};

            double T4 {li8starLV->E() - li8starLV->M()};
            decayStep2->SetDecay(T4, li8starLV->Theta(), li8starLV->Phi());
            double w2 {decayStep2->Generate()};
            nLV = *decayStep2->GetLorentzVector(0);
            alfaLV = *decayStep2->GetLorentzVector(1);
            tritonLV = *decayStep2->GetLorentzVector(2);

            weight = w1 * w2;
        }
        else if(!kUseResonance8Li && kUseResonance7Li)
        {
            kinGen->SetBeamAndExEnergies(TBeam, Ex7Li);
            double w1 {kinGen->Generate()};
            pLV = *kinGen->GetLorentzVector(0);
            auto* li7starLV {kinGen->GetLorentzVector(1)};
            nLV = *kinGen->GetLorentzVector(2); // ya sale directo, sin DecayGenerator

            double T7 {li7starLV->E() - li7starLV->M()};
            decayStep2->SetDecay(T7, li7starLV->Theta(), li7starLV->Phi());
            double w2 {decayStep2->Generate()};
            alfaLV = *decayStep2->GetLorentzVector(0);
            tritonLV = *decayStep2->GetLorentzVector(1);

            weight = w1 * w2;
        }
        else // !kUseResonance8Li && !kUseResonance7Li
        {
            weight = directGen4Body->Generate();
            pLV = *directGen4Body->GetDecay(0);
            nLV = *directGen4Body->GetDecay(1);
            alfaLV = *directGen4Body->GetDecay(2);
            tritonLV = *directGen4Body->GetDecay(3);
            // TGenPhaseSpace trabaja en GeV con las masas que le dimos; reescalar a MeV
            pLV *= 1. / ActPhysics::Constants::kMeVToGeV;
            nLV *= 1. / ActPhysics::Constants::kMeVToGeV;
            alfaLV *= 1. / ActPhysics::Constants::kMeVToGeV;
            tritonLV *= 1. / ActPhysics::Constants::kMeVToGeV;
        }

        // Fill genericos
        h["p"]->Fill(pLV, weight);
        h["n"]->Fill(nLV, weight);
        h["alfa"]->Fill(alfaLV, weight);
        h["t"]->Fill(tritonLV, weight);

        // Balance de energia total (chequeo de validez, valido en los 4 modos)
        double Efinal {pLV.E() + nLV.E() + alfaLV.E() + tritonLV.E()};
        hBalance->Fill(Einitial - Efinal, weight);
    }

    // ---- Dibujar ----
    auto* cTheta = new TCanvas("cTheta", "Angulos Lab", 1200, 800);
    cTheta->Divide(2, 2);
    int pad {1};
    for(const auto& name : {"p", "n", "alfa", "t"})
    {
        cTheta->cd(pad++);
        h[name]->theta->DrawClone("hist");
    }

    auto* cEnergy = new TCanvas("cEnergy", "Energias Lab", 1200, 800);
    cEnergy->Divide(2, 2);
    pad = 1;
    for(const auto& name : {"p", "n", "alfa", "t"})
    {
        cEnergy->cd(pad++);
        h[name]->energy->DrawClone("hist");
    }

    auto* cKin = new TCanvas("cKin", "Cinematica Lab", 1200, 800);
    cKin->Divide(2, 2);
    pad = 1;
    for(const auto& name : {"p", "n", "alfa", "t"})
    {
        cKin->cd(pad++);
        h[name]->kin->DrawClone("colz");
    }

    auto* cBalance = new TCanvas("cBalance", "Balance de energia", 800, 600);
    hBalance->DrawClone("hist");

    std::cout << BOLDGREEN << "Modo simulado: kUseResonance8Li=" << kUseResonance8Li
              << ", kUseResonance7Li=" << kUseResonance7Li << RESET << '\n';
}