#ifndef PIPE3_FILTER_H
#define PIPE3_FILTER_H
#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"
#include "TCanvas.h"
#include "TH1.h"
#include "TMath.h"
#include "TString.h"

#include <fstream>
#include <string>
#include <iostream>

#include "../HistConfig.h"

// ======================================================
// Función auxiliar opcional: guarda los eventos filtrados
// ======================================================
void WriteRejectedEvents(const std::string& infile)
{
    ROOT::DisableImplicitMT(); // Seguridad al escribir

    ROOT::RDataFrame df("Final_Tree", infile);

    // Filtro 1: carga promedio muy alta
    auto rejectedHighQave = df.Filter(
        [](ActRoot::MergerData& m) { return m.fHeavy.fQave > 2000.; },
        {"MergerData"});

    // Filtro 2: no llega al final de ACTAR
    auto rejectedShort = df.Filter(
        [](ActRoot::MergerData& m)
        {
            auto TL = m.fHeavy.fTL;
            auto theta = m.fThetaHeavy * TMath::DegToRad();
            auto z_end = m.fRP.X() + TL * TMath::Cos(theta);
            return z_end < 240.;
        },
        {"MergerData"});

    // Filtro 3: voxeles con carga excesiva
    auto rejectedHighVox = df.Filter(
        [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
        {
            double qTooHigh = 4000.;
            auto heavyCluster = tpc.fClusters[m.fHeavyIdx];
            auto voxelesHeavy = heavyCluster.GetRefToVoxels();
            int counter = 0;
            for (auto& voxel : voxelesHeavy)
                if (voxel.GetCharge() > qTooHigh)
                    counter++;
            return counter > 5;
        },
        {"MergerData", "TPCData"});

    // --- Escribir resultados ---
    std::ofstream out1("./Outputs/rejected_highQave.dat");
    rejectedHighQave.Foreach([&](ActRoot::MergerData& m){ m.Stream(out1); }, {"MergerData"});
    out1.close();

    std::ofstream out2("./Outputs/rejected_shortTrack.dat");
    rejectedShort.Foreach([&](ActRoot::MergerData& m){ m.Stream(out2); }, {"MergerData"});
    out2.close();

    std::ofstream out3("./Outputs/rejected_highVoxelCharge.dat");
    rejectedHighVox.Foreach([&](ActRoot::MergerData& m){ m.Stream(out3); }, {"MergerData"});
    out3.close();
}


void Pipe3_Filter(const std::string& beam, const std::string& target, const std::string& light)
{
    std::string dataconf {};
    if(beam == "11Li")
        dataconf = "./../configs/data.conf";
    else if(beam == "7Li")
        dataconf = "./../configs/data_7Li.conf";
    else
        throw std::runtime_error("Beam cannot differ from 11Li or 7Li");

    // Leer datos
    ActRoot::DataManager dataman {dataconf, ActRoot::ModeType::EReadTPC};

    auto infile {TString::Format("./Outputs/tree_ex_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str())};

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"Final_Tree", infile.Data()};

    // Aplicar filtros
    auto dfFilter = df.Filter(
        [](ActRoot::MergerData& m)
        {
            if(m.fHeavy.fQave > 2000.)
                return false;
            return true;
        },
        {"MergerData"}).Filter(
            [](ActRoot::MergerData& m)
            {
                auto TL {m.fHeavy.fTL};
                auto theta {m.fThetaHeavy * TMath::DegToRad()};
                auto z_end {m.fRP.X() + TL * TMath::Cos(theta)};
                if(z_end < 240)
                    return false;
                return true;
            },
            {"MergerData"}).Filter(
                [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                {
                    double qTooHigh {4000.};
                    auto heavyCluster {tpc.fClusters[m.fHeavyIdx]};
                    auto voxelesHeavy {heavyCluster.GetRefToVoxels()};
                    int counter {};
                    for(auto& voxel : voxelesHeavy)
                    {
                        if(voxel.GetCharge() > qTooHigh)
                            counter++;
                    }
                    if(counter > 5)
                        return false;
                    return true;
                },
                {"MergerData", "TPCData"});

    // Guardar resultado final
    auto outfile {TString::Format("./Outputs/tree_ex_%s_%s_%s_filtered.root", beam.c_str(), target.c_str(), light.c_str())};
    dfFilter.Snapshot("Final_Tree", outfile);
    std::cout << "Saving Final_Tree in " << outfile << '\n';

    // Comparar histogramas
    auto hExBefore = df.Histo1D({"hExBefore", "Excitation Energy before filtering;Ex (MeV);Counts", 100, -5, 10}, "Ex");
    auto hExAfter = dfFilter.Histo1D({"hExAfter", "Excitation Energy after filtering;Ex (MeV);Counts", 100, -5, 10}, "Ex");

    auto c = new TCanvas("cExFilter", "cExFilter", 800, 600);
    c->Divide(2,1);
    c->cd(1);
    hExBefore->DrawClone();
    c->cd(2);
    hExAfter->DrawClone();

    // Opcional: guardar eventos rechazados
    WriteRejectedEvents(infile.Data());
}

#endif
