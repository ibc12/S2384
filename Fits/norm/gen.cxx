#include "PhysExperiment.h"

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

namespace
{
    // Lee un fichero de dos columnas (run  valor) en un map<run, valor>
    std::map<int, double> ReadRunValueFile(const std::string& filename)
    {
        std::ifstream streamer {filename};
        if(!streamer)
            throw std::runtime_error("No se pudo abrir el fichero " + filename);

        std::map<int, double> data {};
        std::string line {};
        while(std::getline(streamer, line))
        {
            if(line.empty())
                continue;
            std::istringstream iss {line};
            int run {};
            double value {};
            if(iss >> run >> value)
                data[run] = value;
        }
        return data;
    }

    // Suma counts/ratio de todos los runs de un beam y genera UN unico
    // fichero de normalizacion con el total corregido
    void GenTotalCorrectedNorm(const std::string& beam, double Nd, double Ndiv,
                                const std::map<int, double>& ratioPerRun)
    {
        auto countsPerRun {ReadRunValueFile("./cfa_perRun_" + beam + ".dat")};

        double totalCorrectedCounts {0};
        for(const auto& [run, counts] : countsPerRun)
        {
            auto ratioIt {ratioPerRun.find(run)};
            if(ratioIt == ratioPerRun.end())
            {
                std::cout << "-> Aviso: run " << run << " no tiene ratio calculado, se omite" << '\n';
                continue;
            }
            double ratio {ratioIt->second};
            if(ratio == 0)
            {
                std::cout << "-> Aviso: run " << run << " tiene ratio 0, se omite" << '\n';
                continue;
            }
            totalCorrectedCounts += counts / ratio;
        }

        PhysUtils::Experiment normBeam {Nd, totalCorrectedCounts, Ndiv};
        std::string filename {"./" + beam + "_norm_L1.dat"};
        normBeam.Print();
        normBeam.Write(filename);
    }
} // namespace

// Generates the different normalizations
void gen()
{
    // Number of incoming beams
    double Ntrigger_7Li {645235};
    double Ntrigger_11Li {3079721};
    double Ndiv {300};
    // Total ACTAR length used in the LISE calculation
    double totalLength {25.6}; // cm
    double Nd {4.126e19 * totalLength}; // number of deuterons in the target

    // Build
    PhysUtils::Experiment norm11Li {Nd, Ntrigger_11Li, Ndiv};
    norm11Li.Print();
    norm11Li.Write("./11Li_norm.dat");

    PhysUtils::Experiment norm7Li {Nd, Ntrigger_7Li, Ndiv};
    norm7Li.Print();
    norm7Li.Write("./7Li_norm.dat");

    // El fichero de ratios CFA/F2 es comun a ambos haces: se lee una unica vez
    auto ratioPerRun {ReadRunValueFile("./CFA_F2_ratio_run.dat")};

    // Normalizacion total corregida por beam, sumando counts/ratio de cada run
    GenTotalCorrectedNorm("11Li", Nd, Ndiv, ratioPerRun);
    GenTotalCorrectedNorm("7Li", Nd, Ndiv, ratioPerRun);
}