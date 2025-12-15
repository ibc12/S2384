#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActMergerData.h"

#include "ROOT/RDataFrame.hxx"

#include <fstream>
#include <string>

void EventsSelectionPIDFront()
{
    // ------------------------
    // Configuración de datos
    // ------------------------
    std::string beam {"11Li"};
    std::string dataconf {};
    if(beam == "11Li")
        dataconf = "./../configs/data_11Li.conf";
    else if(beam == "7Li")
        dataconf = "./../configs/data_7Li.conf";
    else
        throw std::runtime_error("Beam must be 11Li or 7Li");

    ActRoot::DataManager dataman {dataconf, ActRoot::ModeType::EMerge};

    auto chain = dataman.GetChain();
    auto chain2 = dataman.GetChain(ActRoot::ModeType::EReadSilMod);
    auto chain3 = dataman.GetChain(ActRoot::ModeType::EFilter);
    auto chain4 = dataman.GetChain(ActRoot::ModeType::EReadTPC);

    chain->AddFriend(chain2.get());
    chain->AddFriend(chain3.get(), "TPCData");
    chain->AddFriend(chain4.get(), "GETTree");

    //ROOT::EnableImplicitMT();
    ROOT::RDataFrame df0 {*chain};

    // -----------------------------------
    // Filtro de pads / runs (IGUAL QUE PIPE1)
    // -----------------------------------
    auto dfPads = df0.Filter(
        [](ActRoot::MergerData& m)
        {
            if(m.fLight.fNs.size())
            {
                if(m.fLight.fLayers.front() == "l0" && m.fLight.fNs.front() == 9)
                    return false;
                if(m.fLight.fLayers.front() == "f0" && m.fLight.fNs.front() == 2)
                    return false;
                if(m.fLight.fLayers.front() == "r0" && m.fLight.fNs.front() == 3)
                    return false;
            }

            if(m.fRun > 29 && m.fRun < 35)
            {
                if(m.fRun == 30 || m.fRun == 31 || m.fRun == 33)
                {
                    if(!m.fLight.fLayers.empty() && m.fLight.fLayers.front() == "r0")
                    {
                        if(!m.fLight.fNs.empty() && (m.fLight.fNs.front() == 3 || m.fLight.fNs.front() == 5))
                            return false;
                    }
                }
                if(m.fRun == 32)
                {
                    if(!m.fLight.fLayers.empty() && m.fLight.fLayers.front() == "r0")
                        return false;
                }
                if(m.fRun == 34)
                {
                    if(!m.fLight.fLayers.empty() &&
                       (m.fLight.fLayers.front() == "r0" || m.fLight.fLayers.front() == "l0"))
                        return false;
                }
            }

            if(m.fRun > 35 && m.fRun < 45)
            {
                if(!m.fLight.fLayers.empty() && m.fLight.fLayers.front() == "f0" && !m.fLight.fNs.empty() &&
                   m.fLight.fNs.front() == 5)
                    return false;
            }

            if(m.fRun == 116 || m.fRun == 117)
            {
                if(!m.fLight.fLayers.empty() && m.fLight.fLayers.front() == "r0" && !m.fLight.fNs.empty() &&
                   m.fLight.fNs.front() == 2)
                    return false;
            }

            return true;
        },
        {"MergerData"});

    // ------------------------
    // Evento válido:
    //  - 1 solo silicio
    //  - capa f0
    // ------------------------
    auto isOneSilF0 = [](ActRoot::MergerData& m)
    { return (m.fLight.GetNLayers() == 1 && m.fLight.GetLayer(0) == "f0"); };

    // ------------------------
    // Leer corte PID f0
    // ------------------------
    ActRoot::CutsManager<std::string> cuts;
    cuts.ReadCut("f0", TString::Format("./Cuts/events3He_11Li.root").Data());

    // ------------------------
    // Aplicar PID
    // ------------------------

    auto gated = dfPads.Filter(
        [&](ActRoot::MergerData& m)
        {
            if(!isOneSilF0(m))
                return false;

            if(!cuts.GetCut("f0"))
                return false;

            return cuts.IsInside("f0", m.fLight.fEs[0], m.fLight.fQave);
        },
        {"MergerData"});

    // ------------------------
    // Guardar eventos
    // ------------------------
    std::ofstream out(TString::Format("./Outputs/events_3He_f0.dat").Data());

    gated.Foreach([&](ActRoot::MergerData& m) { m.Stream(out); }, {"MergerData"});

    out.close();

    std::cout << "Eventos f0 guardados correctamente." << std::endl;
}
