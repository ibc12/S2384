// ReconstructExElasticAsProtonL1.cxx
//
// Lee el tree de Ex, se queda con los eventos con Ex en (-1, 1) MeV
// (canal elastico, deuterones) y reconstruye Ex y ThetaCM tratando la particula
// ligera como si fuese un proton:
//   - Energia en el vertice con la tabla de SRIM del proton (1H)
//   - Cinematica 7Li(d,p)8Li
//
// Nuevas columnas: Ex_d_as_p y ThetaCM_d_as_p
//
// Uso: root -l -b -q ReconstructExElasticAsProtonL1.cxx

#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActParticle.h"
#include "ActSRIM.h"

#include "ROOT/RDataFrame.hxx"

#include "TMath.h"
#include "TString.h"

#include <iostream>
#include <string>
#include <vector>

void ReconstructExElasticAsProtonL1()
{
    // ---------------- Configuracion ----------------
    const std::string inFile {"../../PostAnalysis/Outputs/tree_ex_F_7Li_d_d.root"};
    const std::string outFile {"./Outputs/tree_ex_F_7Li_d_as_p.root"};
    const std::string treeName {"Final_Tree"}; // <-- cambiar si el TTree tiene otro nombre

    const std::string beam {"7Li"};
    const std::string target {"2H"};
    const std::string light {"p"}; // la particula se trata COMO proton

    const double exMin {-1.}; // MeV
    const double exMax {1.};  // MeV

    // ---------------- SRIM ----------------
    auto* srim {new ActPhysics::SRIM};
    // Tabla del proton: se usa para la energia en el vertice
    srim->ReadTable(light, "../../Calibrations/SRIM/1H_900mb_CF4_95-5.txt");
    srim->ReadTable(beam, TString::Format("../../Calibrations/SRIM/%s_900mb_CF4_95-5.txt", beam.c_str()).Data());
    srim->ReadTable("mylar", TString::Format("../../Calibrations/SRIM/%s_Mylar.txt", beam.c_str()).Data());

    // ---------------- Particulas ----------------
    ActPhysics::Particle pb {beam};
    ActPhysics::Particle pt {target};
    ActPhysics::Particle pl {light};

    // ---------------- Energia inicial del haz ----------------
    double initialEnergy {7.558}; // medida por los operadores; resolucion del 0.19%
    initialEnergy = srim->Slow("mylar", initialEnergy * pb.GetAMU(), 0.0168);
    initialEnergy = srim->Slow(beam, initialEnergy, 60); // 60 mm de gas antes del pad plane
    initialEnergy = initialEnergy / pb.GetAMU();         // de vuelta a MeV/u
    std::cout << "Using initial energy for Pipe2: " << initialEnergy << " AMeV" << std::endl;

    // Cinematica 7Li(d,p)8Li
    ActPhysics::Kinematics kin {pb, pt, pl, initialEnergy * pb.GetAMU()};

    // ---------------- DataFrame ----------------
    ROOT::RDataFrame df {treeName, inFile};

    // Solo el intervalo Ex in (-1, 1) y solo L1 events (fLight no filled)
    auto dfSel = df.Filter([&](double Ex) { return Ex > exMin && Ex < exMax; }, {"Ex"})
                     .Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"});

    // Energia en el vertice asumiendo que la particula es un proton
    auto dfVertex = dfSel.Define("EVertex_d_as_p",
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

    // Energia del haz en el vertice (no depende de la particula ligera)
    auto dfBeam = dfVertex.Define("EBeam_d_as_p", [&](const ActRoot::MergerData& d)
                                  { return srim->Slow(beam, initialEnergy * pb.GetAMU(), d.fRP.X()); }, {"MergerData"});

    // Un objeto Kinematics por slot (cambiamos EBeam en cada evento)
    std::vector<ActPhysics::Kinematics> vkins {dfBeam.GetNSlots()};
    for(auto& k : vkins)
        k = kin;

    auto dfEx = dfBeam.DefineSlot("Ex_d_as_p",
                                  [&](unsigned int slot, const ActRoot::MergerData& d, double EVertex, double EBeam)
                                  {
                                      vkins[slot].SetBeamEnergy(EBeam);
                                      return vkins[slot].ReconstructExcitationEnergy(EVertex, (d.fThetaLight) *
                                                                                                  TMath::DegToRad());
                                  },
                                  {"MergerData", "EVertex_d_as_p", "EBeam_d_as_p"});

    auto dfFinal = dfEx.DefineSlot("ThetaCM_d_as_p",
                                   [&](unsigned int slot, const ActRoot::MergerData& d, double EVertex, double EBeam)
                                   {
                                       vkins[slot].SetBeamEnergy(EBeam);
                                       return vkins[slot].ReconstructTheta3CMFromLab(EVertex, (d.fThetaLight) *
                                                                                                  TMath::DegToRad()) *
                                              TMath::RadToDeg();
                                   },
                                   {"MergerData", "EVertex_d_as_p", "EBeam_d_as_p"});

    // ---------------- Guardar ----------------
    // Se guardan todas las columnas originales + las nuevas
    dfFinal.Snapshot(treeName, outFile);
    std::cout << "Tree guardado en: " << outFile << std::endl;
}