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

    struct BeamSums
    {
        double raw {0};       // suma de counts (sin corregir)
        double corrected {0}; // suma de counts / ratio
        int nRuns {0};        // runs del intervalo
    };

    // Suma los runs de un beam que caen dentro del intervalo
    BeamSums SumRuns(const std::string& beam, const std::map<int, double>& ratioPerRun,
                     int runMin, int runMax)
    {
        auto countsPerRun {ReadRunValueFile("../cfa_perRun_" + beam + ".dat")};
        BeamSums sums {};
        for(const auto& [run, counts] : countsPerRun)
        {
            if(!InRunRange(run, runMin, runMax))
                continue;
            sums.raw += counts;
            sums.nRuns++;

            auto ratioIt {ratioPerRun.find(run)};
            if(ratioIt == ratioPerRun.end())
            {
                std::cout << "-> Aviso: run " << run << " no tiene ratio calculado, se omite en la corregida" << '\n';
                continue;
            }
            if(ratioIt->second == 0)
            {
                std::cout << "-> Aviso: run " << run << " tiene ratio 0, se omite en la corregida" << '\n';
                continue;
            }
            sums.corrected += counts / ratioIt->second;
        }
        return sums;
    }

    // Genera las dos normalizaciones de un beam: sin corregir y corregida
    void GenNorms(const std::string& beam, double Ntrigger, double Nd, double Ndiv,
                  const std::map<int, double>& ratioPerRun, int runMin, int runMax)
    {
        auto sums {SumRuns(beam, ratioPerRun, runMin, runMax)};
        bool hasRange {runMin >= 0 || runMax >= 0};
        std::string suffix {RunRangeSuffix(runMin, runMax)};

        if(sums.nRuns == 0)
            std::cout << "-> Aviso: ningun run de " << beam << " dentro del intervalo" << '\n';
        std::cout << "-> " << beam << ": " << sums.nRuns << " runs, suma CFA = " << sums.raw
                  << " (Ntrigger = " << Ntrigger << ")" << '\n';

        // Sin corregir: Ntrigger global si no hay intervalo, suma de CFA si lo hay
        PhysUtils::Experiment norm {Nd, hasRange ? sums.raw : Ntrigger, Ndiv};
        norm.Print();
        norm.Write("./" + beam + "_norm" + suffix + ".dat");

        // Corregida
        PhysUtils::Experiment normCorr {Nd, sums.corrected, Ndiv};
        normCorr.Print();
        normCorr.Write("./" + beam + "_norm_L1" + suffix + ".dat");
    }
} // namespace

// Generates the different normalizations
// runMin / runMax: intervalo (inclusivo) de runs. Con -1 el limite no se aplica.
void gen(int runMin = 59, int runMax = 122)
{
    // Number of incoming beams
    double Ntrigger_7Li {645235};
    double Ntrigger_11Li {3079721};
    double Ndiv {300};
    // Total ACTAR length used in the LISE calculation
    double totalLength {25.6}; // cm
    double Nd {4.126e19 * totalLength}; // number of deuterons in the target

    // El fichero de ratios CFA/F2 es comun a ambos haces: se lee una unica vez
    auto ratioPerRun {ReadRunValueFile("../CFA_F2_ratio_run.dat")};

    GenNorms("11Li", Ntrigger_11Li, Nd, Ndiv, ratioPerRun, runMin, runMax);
    GenNorms("7Li", Ntrigger_7Li, Nd, Ndiv, ratioPerRun, runMin, runMax);
}