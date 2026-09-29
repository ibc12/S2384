#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActParticle.h"
#include "ActSRIM.h"
#include "ActSilData.h"
#include "ActSilMatrix.h"
#include "ActSilSpecs.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TF1.h"
#include "TH2.h"
#include "TH2D.h"
#include "TPaveText.h"
#include "TString.h"

#include <algorithm>
#include <fstream>
#include <map>
#include <string>
#include <utility>

#include "../../../Fits//Histos.h"
#include "../../../PostAnalysis/HistConfig.h"
#include "../../../PrettyStyle.C"
#include "../Utils.h"


// This macro is to check that multiplicity 3 evetns fall on the alfa+t PID for the front silicon

struct DecayInfo
{
    int maxThetaIdx;
    int minThetaIdx;
    int maxQLengthIdx;
    int minQLengthIdx;

    double maxTheta;
    double minTheta;
    double maxQLength;
    double minQLength;

    double beamQLength;
    double lightQLength;
};

// Resultado del matching geometrico particula<->silicio, ya organizado en los dos telescopios frontales.
// Cada par es {E_capa_fina, E_capa_gruesa}, es decir {f0,f1} o {f2,f3}.
// Convencion: {0,0} significa "sin choque valido en este telescopio para esta particula"
// (bien porque no hubo ningun hit, bien porque solo hubo hit en la capa gruesa sin la fina).
struct TelescopeSil
{
    std::pair<double, double> alfa01 {0., 0.};   // {E_f0, E_f1} para la alfa
    std::pair<double, double> triton01 {0., 0.}; // {E_f0, E_f1} para el triton
    std::pair<double, double> alfa23 {0., 0.};   // {E_f2, E_f3} para la alfa
    std::pair<double, double> triton23 {0., 0.}; // {E_f2, E_f3} para el triton
    bool alfaConflict {false};                   // true si la alfa matcheo (antes de resolver) en los 2 telescopios
    bool tritonConflict {false};                 // idem para el triton
};


void CheckPID_alfa_triton()
{
    // We had already identify the alfa and triton by the Q/L ratios

    std::string beam {"7Li"};
    std::string target {"d"};
    std::string light {"p"};

    // Get silspecs files for post (7Li is also post)
    auto specs {std::make_shared<ActPhysics::SilSpecs>()};
    std::string silConfig("silspecs");
    specs->ReadFile("../../../configs/" + silConfig + ".conf");

    // Get file from pipe3
    TString infile = TString::Format("../Outputs/DecayM4_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str());
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df("Final_Tree", infile.Data());

    // Read cuts to later plot
    ActRoot::CutsManager<std::string> cuts;
    cuts.ReadCut("d", "../Cuts/pid_d_f0_11Li.root");
    cuts.ReadCut("t", "../Cuts/pid_t_f0_11Li.root");
    cuts.ReadCut("3He", "../Cuts/pid_3He_f0_11Li.root");
    cuts.ReadCut("a", "../Cuts/pid_4He_f0_11Li.root");

    // Define alfaidx and triton idx variables
    auto def = df.Define("AlfaIdx", [](DecayInfo decay) { return decay.maxQLengthIdx; }, {"Decay"})
                   .Define("TritonIdx", [](DecayInfo decay) { return decay.minQLengthIdx; }, {"Decay"});

    // --------------------------------------------------------------------------------------
    // Helpers geometricos para hacer el matching entre las trazas del TPC y los hits de silicio
    // --------------------------------------------------------------------------------------

    // Tolerancia maxima de matching (mismo criterio que en Pipe1_PIDM4, kSilMatchTol)
    const double kMatchTol {170.}; // mm

    // Layers del wall frontal, agrupadas por telescopio: (f0,f1) primer telescopio, (f2,f3) segundo.
    // Se asume que la primera layer de cada pareja (f0, f2) es la capa fina/dE de entrada.
    const std::vector<std::string> frontLayers {"f0", "f1", "f2", "f3"};

    // Punto de impacto esperado (proyeccion) de la traza de una particula sobre el plano de una layer,
    // siguiendo el mismo procedimiento que Pipe1_PIDM4 usa para "f0"
    auto getProjection = [&](ActRoot::TPCData& tpc, int idx, const std::string& layer) -> ROOT::Math::XYZPointF
    {
        auto& cluster = tpc.fClusters[idx];
        auto line = cluster.GetLine();
        line.Scale(Utils::scaleXY, Utils::scaleZ);
        auto pointLayer = specs->GetLayer(layer).GetPoint();
        return line.MoveToX(pointLayer.X());
    };

    // Posicion real de un pad concreto dentro de una layer.
    // NOTA: se asume la misma convencion que usa Pipe1_PIDM4 para "f0": el punto de la layer
    // fija la profundidad (X) y GetPlacements() da la posicion (Y,Z) del pad en ese plano.
    // Ajustar aqui si la convencion real de SilSpecs/SilLayer es distinta.
    auto getHitPos = [&](const std::string& layer, int padIdx) -> ROOT::Math::XYZPoint
    {
        auto pointLayer = specs->GetLayer(layer).GetPoint();
        auto placement = specs->GetLayer(layer).GetPlacements().at(padIdx);
        return ROOT::Math::XYZPoint(pointLayer.X(), placement.first, placement.second);
    };

    // Construye el par {dE, E} de un telescopio aplicando la regla de validez:
    // solo cuenta si hay hit en la capa fina (dE); si no, el evento no es valido para ese telescopio.
    auto buildTelescopePair = [](double dE, double E) -> std::pair<double, double>
    {
        if(dE > 0.)
            return {dE, (E > 0. ? E : 0.)};
        return {0., 0.};
    };

    // Now get if there is any hit of the silicon corresponding to any of the particles
    // Diferent cases - 2 front; 1 front; 1 telescope; nothing
    // Let's treat first the case of the front layers and then the telescope layers
    auto defTel =
        def.Define("TelInfo",
                   [&](int alfaIdx, int tritonIdx, ActRoot::TPCData& tpc, ActRoot::SilData& sil) -> TelescopeSil
                   {
                       TelescopeSil res {};

                       if(alfaIdx < 0 || tritonIdx < 0)
                           return res;

                       // First check if there is silicon hit in any of the layers applying the threshold
                       sil.ApplyFinerThresholds(specs);

                       auto layersHit = sil.GetLayers();
                       bool hasFront = std::any_of(layersHit.begin(), layersHit.end(), [](const std::string& l)
                                                   { return l == "f0" || (l == "f2" && l == "f3"); });
                       if(!hasFront)
                           return res;

                       // If it has front check if more than one impact:
                       // proyeccion esperada de cada particula sobre cada una de las 4 layers frontales
                       std::map<std::string, ROOT::Math::XYZPoint> alfaProj, tritonProj;
                       for(const auto& layer : frontLayers)
                       {
                           alfaProj[layer] = getProjection(tpc, alfaIdx, layer);
                           tritonProj[layer] = getProjection(tpc, tritonIdx, layer);
                       }

                       // Energia y distancia de matching asignadas a cada particula, por layer.
                       // -1 en energia / distancia enorme => no se le asigno ningun hit en esa layer.
                       std::map<std::string, double> alfaE, tritonE, alfaDist, tritonDist;
                       for(const auto& layer : frontLayers)
                       {
                           alfaE[layer] = -1.;
                           tritonE[layer] = -1.;
                           alfaDist[layer] = 1e9;
                           tritonDist[layer] = 1e9;
                       }

                       for(const auto& layerStr : frontLayers)
                       {
                           auto it = sil.fSiN.find(layerStr);
                           if(it == sil.fSiN.end() || it->second.empty())
                               continue; // no hay hit en esta layer

                           const auto& pads = it->second;                // indices de los pads que dispararon
                           const auto& energies = sil.fSiE.at(layerStr); // energias, mismo orden que pads

                           if(pads.size() == 1)
                           {
                               // Solo un hit: puede ser de la alfa o del triton (para f2/f3 solo deberia llegar
                               // una de las dos, pero igualmente decidimos por distancia)
                               auto hitPos = getHitPos(layerStr, pads[0]);
                               double dAlfa = (hitPos - alfaProj[layerStr]).R();
                               double dTriton = (hitPos - tritonProj[layerStr]).R();

                               if(dAlfa <= dTriton && dAlfa < kMatchTol)
                               {
                                   alfaE[layerStr] = energies[0];
                                   alfaDist[layerStr] = dAlfa;
                               }
                               else if(dTriton < dAlfa && dTriton < kMatchTol)
                               {
                                   tritonE[layerStr] = energies[0];
                                   tritonDist[layerStr] = dTriton;
                               }
                           }
                           else
                           {
                               // Dos (o mas) hits: nos quedamos con los dos primeros y resolvemos como una
                               // asignacion 2 a 2, probando las dos combinaciones posibles y quedandonos con
                               // la de menor distancia total
                               int i0 {0}, i1 {1};
                               auto hitPos0 = getHitPos(layerStr, pads[i0]);
                               auto hitPos1 = getHitPos(layerStr, pads[i1]);

                               double dAlfa0 = (hitPos0 - alfaProj[layerStr]).R();
                               double dAlfa1 = (hitPos1 - alfaProj[layerStr]).R();
                               double dTriton0 = (hitPos0 - tritonProj[layerStr]).R();
                               double dTriton1 = (hitPos1 - tritonProj[layerStr]).R();

                               double totalCombo1 = dAlfa0 + dTriton1; // alfa->hit0, triton->hit1
                               double totalCombo2 = dAlfa1 + dTriton0; // alfa->hit1, triton->hit0

                               if(totalCombo1 <= totalCombo2)
                               {
                                   if(dAlfa0 < kMatchTol)
                                   {
                                       alfaE[layerStr] = energies[i0];
                                       alfaDist[layerStr] = dAlfa0;
                                   }
                                   if(dTriton1 < kMatchTol)
                                   {
                                       tritonE[layerStr] = energies[i1];
                                       tritonDist[layerStr] = dTriton1;
                                   }
                               }
                               else
                               {
                                   if(dAlfa1 < kMatchTol)
                                   {
                                       alfaE[layerStr] = energies[i1];
                                       alfaDist[layerStr] = dAlfa1;
                                   }
                                   if(dTriton0 < kMatchTol)
                                   {
                                       tritonE[layerStr] = energies[i0];
                                       tritonDist[layerStr] = dTriton0;
                                   }
                               }
                           }
                       }

                       // Pares dE-E por telescopio y por particula
                       res.alfa01 = buildTelescopePair(alfaE["f0"], alfaE["f1"]);
                       res.triton01 = buildTelescopePair(tritonE["f0"], tritonE["f1"]);
                       res.alfa23 = buildTelescopePair(alfaE["f2"], alfaE["f3"]);
                       res.triton23 = buildTelescopePair(tritonE["f2"], tritonE["f3"]);

                       // Comprobacion: una misma particula no deberia chocar en los dos telescopios a la vez.
                       // Si ocurre (por ejemplo por un mal matching), nos quedamos con el telescopio cuya capa
                       // fina tuvo mejor (menor) distancia de matching, y anulamos el otro.
                       res.alfaConflict = (res.alfa01.first > 0.) && (res.alfa23.first > 0.);
                       res.tritonConflict = (res.triton01.first > 0.) && (res.triton23.first > 0.);

                       if(res.alfaConflict)
                       {
                           if(alfaDist["f0"] <= alfaDist["f2"])
                               res.alfa23 = {0., 0.};
                           else
                               res.alfa01 = {0., 0.};
                       }
                       if(res.tritonConflict)
                       {
                           if(tritonDist["f0"] <= tritonDist["f2"])
                               res.triton23 = {0., 0.};
                           else
                               res.triton01 = {0., 0.};
                       }

                       return res;
                   },
                   {"AlfaIdx", "TritonIdx", "TPCData", "SilData"});

    // Columnas individuales, mas comodas de usar/plotear
    auto defPairs =
        defTel.Define("alfaTelescope01", [](const TelescopeSil& t) { return t.alfa01; }, {"TelInfo"})
            .Define("tritonTelescope01", [](const TelescopeSil& t) { return t.triton01; }, {"TelInfo"})
            .Define("alfaTelescope", [](const TelescopeSil& t) { return t.alfa23; }, {"TelInfo"})
            .Define("tritonTelescope", [](const TelescopeSil& t) { return t.triton23; }, {"TelInfo"})
            .Define("AlfaTelescopeConflict", [](const TelescopeSil& t) { return t.alfaConflict; }, {"TelInfo"})
            .Define("TritonTelescopeConflict", [](const TelescopeSil& t) { return t.tritonConflict; }, {"TelInfo"});

    // Energia total por particula, derivada del telescopio (unico) que resulto valido tras la
    // comprobacion anterior. -1 si no hubo ningun choque valido.
    auto defE = defPairs
                    .Define("AlfaE",
                            [](const std::pair<double, double>& t01, const std::pair<double, double>& t23)
                            {
                                if(t01.first > 0.)
                                    return t01.first + t01.second;
                                if(t23.first > 0.)
                                    return t23.first + t23.second;
                                return -1.;
                            },
                            {"alfaTelescope01", "alfaTelescope"})
                    .Define("TritonE",
                            [](const std::pair<double, double>& t01, const std::pair<double, double>& t23)
                            {
                                if(t01.first > 0.)
                                    return t01.first + t01.second;
                                if(t23.first > 0.)
                                    return t23.first + t23.second;
                                return -1.;
                            },
                            {"tritonTelescope01", "tritonTelescope"});

    // Chequeo rapido: PID de energias de silicio frontal para eventos con ambas particulas asignadas
    auto dfBoth = defE.Filter([](double a, double t) { return a > 0 && t > 0; }, {"AlfaE", "TritonE"});

    // PID plot Qave/Qlength vs SilE de f0, por separado para cada particula (equivalente a hPIDf0
    // de Pipe1_PIDM4, pero usando solo la energia de la capa fina f0 de cada telescopio y el
    // QLength correspondiente a cada particula, guardado en el Decay original)
    auto defQL =
        defE.Define("AlfaF0E", [](const std::pair<double, double>& p) { return p.first; }, {"alfaTelescope01"})
            .Define("TritonF0E", [](const std::pair<double, double>& p) { return p.first; }, {"tritonTelescope01"})
            .Define("AlfaF2E", [](const std::pair<double, double>& p) { return p.first; }, {"alfaTelescope01"})
            .Define("TritonF2E", [](const std::pair<double, double>& p) { return p.first; }, {"tritonTelescope01"})
            .Define("AlfaF3E", [](const std::pair<double, double>& p) { return p.second; }, {"alfaTelescope"})
            .Define("TritonF3E", [](const std::pair<double, double>& p) { return p.second; }, {"tritonTelescope"})
            .Define("AlfaQLength", [](DecayInfo decay) { return decay.maxQLength; }, {"Decay"})
            .Define("TritonQLength", [](DecayInfo decay) { return decay.minQLength; }, {"Decay"});

    // Solo tiene sentido pintar el punto si la particula realmente choco en f0
    auto dfAlfaF0 = defQL.Filter([](double e) { return e > 0; }, {"AlfaF0E"});
    auto dfTritonF0 = defQL.Filter([](double e) { return e > 0; }, {"TritonF0E"});

    // NOTA: el rango del eje Y (QLength) es una estimacion; ajustalo a los valores tipicos que
    // tengas en tu Decay (minQLength/maxQLength) para que el histograma no salga vacio o saturado
    auto hPIDf0Alfa =
        dfAlfaF0.Histo2D({"hPIDf0Alfa", "PID plot QLength vs SilE (alfa, f0);Silicon Energy f0 (MeV);QLength (a.u.)",
                          100, 0, 60, 100, 0, 3000},
                         "AlfaF0E", "AlfaQLength");
    auto hPIDf0Triton = dfTritonF0.Histo2D(
        {"hPIDf0Triton", "PID plot QLength vs SilE (triton, f0);Silicon Energy f0 (MeV);QLength (a.u.)", 100, 0, 60,
         100, 0, 3000},
        "TritonF0E", "TritonQLength");
    auto hPIDf2f3Alfa = defQL.Histo2D(
        {"hPIDf2f3Alfa", "PID plot telescope f2-f3;Silicon Energy f3 (MeV);#Delta E f2 (MeV)", 100, 0, 12, 100, 0, 60},
        "AlfaF3E", "AlfaF2E");
    auto hPIDf2f3Triton =
        defQL.Histo2D({"hPIDf2f3Triton", "PID plot telescope f2-f3;Silicon Energy f3 (MeV);#Delta E f2 (MeV)", 100, 0,
                       12, 100, 0, 60},
                      "TritonF3E", "TritonF2E");

    TCanvas* cPIDf0 = new TCanvas("cPIDf0", "PID f0: SilE vs QLength", 1600, 600);
    cPIDf0->Divide(2, 1);
    cPIDf0->cd(1);
    hPIDf0Alfa->DrawClone("colz");
    cuts.DrawAll();
    cPIDf0->cd(2);
    hPIDf0Triton->DrawClone("colz");
    cuts.DrawAll();

    auto nConflictAlfa = defPairs.Filter([](bool c) { return c; }, {"AlfaTelescopeConflict"}).Count();
    auto nConflictTriton = defPairs.Filter([](bool c) { return c; }, {"TritonTelescopeConflict"}).Count();
    TCanvas* cPIDf2f3 = new TCanvas("cPIDf2f3", "PID f2+f3: SilE vs QLength", 1600, 600);
    cPIDf2f3->Divide(2, 1);
    cPIDf2f3->cd(1);
    hPIDf2f3Alfa->DrawClone("colz");
    cPIDf2f3->cd(2);
    hPIDf2f3Triton->DrawClone("colz");

    std::cout << "Eventos totales                              : " << defE.Count().GetValue() << '\n';
    std::cout << "Eventos con energia asignada a ambas          : " << dfBoth.Count().GetValue() << '\n';
    std::cout << "Conflictos alfa en los dos telescopios a la vez  : " << nConflictAlfa.GetValue() << '\n';
    std::cout << "Conflictos triton en los dos telescopios a la vez: " << nConflictTriton.GetValue() << '\n';
}