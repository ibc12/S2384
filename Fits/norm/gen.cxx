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

    // Comprueba si un run esta dentro del intervalo (-1 = sin limite)
    bool InRunRange(int run, int runMin, int runMax)
    {
        if(runMin >= 0 && run < runMin)
            return false;
        if(runMax >= 0 && run > runMax)
            return false;
        return true;
    }

    // Sufijo para el nombre del fichero, vacio si no hay intervalo
    std::string RunRangeSuffix(int runMin, int runMax)
    {
        if(runMin < 0 && runMax < 0)
            return "";
        std::string first {runMin >= 0 ? std::to_string(runMin) : "first"};
        std::string last {runMax >= 0 ? std::to_string(runMax) : "last"};
        return "_runs" + first + "-" + last;
    }

    // Suma counts/ratio de los runs de un beam (dentro del intervalo) y genera
    // UN unico fichero de normalizacion con el total corregido
    void GenTotalCorrectedNorm(const std::string& beam, double Nd, double Ndiv,
                               const std::map<int, double>& ratioPerRun,
                               int runMin, int runMax)
    {
        auto countsPerRun {ReadRunValueFile("./cfa_perRun_" + beam + ".dat")};

        double totalCorrectedCounts {0};
        int nUsedRuns {0};
        for(const auto& [run, counts] : countsPerRun)
        {
            if(!InRunRange(run, runMin, runMax))
                continue;

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
            nUsedRuns++;
        }

        if(nUsedRuns == 0)
            std::cout << "-> Aviso: ningun run de " << beam << " dentro del intervalo seleccionado" << '\n';
        else
            std::cout << "-> " << beam << ": " << nUsedRuns << " runs usados" << '\n';

        PhysUtils::Experiment normBeam {Nd, totalCorrectedCounts, Ndiv};
        std::string filename {"./" + beam + "_norm_L1" + RunRangeSuffix(runMin, runMax) + ".dat"};
        normBeam.Print();
        normBeam.Write(filename);
    }
} // namespace

// Generates the different normalizations
// runMin / runMax: intervalo (inclusivo) de runs para la normalizacion corregida.
// Con -1 el limite correspondiente no se aplica.
void gen(int runMin = -1, int runMax = -1)
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
    GenTotalCorrectedNorm("11Li", Nd, Ndiv, ratioPerRun, runMin, runMax);
    GenTotalCorrectedNorm("7Li", Nd, Ndiv, ratioPerRun, runMin, runMax);
}