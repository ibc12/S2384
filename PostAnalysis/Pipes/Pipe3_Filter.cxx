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
#include <iostream>
#include <string>

#include "../HistConfig.h"

// ======================================================
// Función auxiliar opcional: guarda los eventos filtrados
// ======================================================
void WriteRejectedEvents(const std::string& infile)
{
    ROOT::DisableImplicitMT(); // Seguridad al escribir

    ROOT::RDataFrame df("Final_Tree", infile);

    // Filtro 1: carga promedio muy alta
    auto rejectedHighQave = df.Filter([](ActRoot::MergerData& m) { return m.fHeavy.fQave > 2000.; }, {"MergerData"});

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

    // Filtro 4: carga en el RP o en el haz
    auto rejectedHighVoxRPBeam = df.Filter(
        [&](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
        {
            auto rp {m.fRP};
            auto rp_y {rp.Y() / 2};
            // Run for all clusters
            int counter {};
            for(auto& cluster : tpc.fClusters)
            {
                auto voxels {cluster.GetRefToVoxels()};
                for(auto& v : voxels)
                {
                    if(v.GetPosition().Y() > rp_y - 3 && v.GetPosition().Y() < rp_y + 3) // aprox L1 exclusion zone
                        if(v.GetCharge() > 3000.)
                            counter++;
                }
            }
            return counter > 4;
        },
        {"MergerData", "GETTree_TPCData"});

    // --- Escribir resultados ---
    std::ofstream out1("./Outputs/rejected_highQave.dat");
    rejectedHighQave.Foreach([&](ActRoot::MergerData& m) { m.Stream(out1); }, {"MergerData"});
    out1.close();

    std::ofstream out2("./Outputs/rejected_shortTrack.dat");
    rejectedShort.Foreach([&](ActRoot::MergerData& m) { m.Stream(out2); }, {"MergerData"});
    out2.close();

    std::ofstream out4("./Outputs/rejected_highVoxelRPBeam.dat");
    rejectedHighVoxRPBeam.Foreach([&](ActRoot::MergerData& m) { m.Stream(out4); }, {"MergerData"});
    out4.close();
}


void Pipe3_Filter(const std::string& beam, const std::string& target, const std::string& light)
{
    bool onlySil {false};

    auto infile {TString::Format("./Outputs/tree_ex_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str())};

    // ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"Final_Tree", infile.Data()};

    // Aplicar filtros
    auto dfFilter = df.Filter( // Check for heavier clusters than Li
                          [](ActRoot::MergerData& m)
                          {
                              if(m.fHeavy.fQave > 2000.)
                                  return false;
                              return true;
                          },
                          {"MergerData"})
                        .Filter( // Check if heavy reaches end of ACTAR (other way to mask heavier clusters,
                                 // but can delete good events with heavy particle bad reconstructed)
                            [](ActRoot::MergerData& m)
                            {
                                auto TL {m.fHeavy.fTL};
                                auto theta {m.fThetaHeavy * TMath::DegToRad()};
                                auto z_end {m.fRP.X() + TL * TMath::Cos(theta)};
                                if(z_end < 240)
                                    return false;
                                return true;
                            },
                            {"MergerData"})
                        .Filter( // Most of the times high charge deposit
                                 // is masked by rp or is in beam cluster
                            [&](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                            {
                                auto rp {m.fRP};
                                auto rp_y {rp.Y() / 2};
                                // Run for all clusters
                                int counter {};
                                for(auto& cluster : tpc.fClusters)
                                {
                                    auto voxels {cluster.GetRefToVoxels()};
                                    for(auto& v : voxels)
                                    {
                                        if(v.GetPosition().Y() > rp_y - 3 &&
                                           v.GetPosition().Y() < rp_y + 3) // aprox L1 exclusion zone
                                            if(v.GetCharge() > 3000.)
                                                counter++;
                                    }
                                }
                                if(counter > 4)
                                    return false;
                                return true;
                            },
                            {"MergerData", "GETTree_TPCData"});

    // Guardar resultado final
    auto outfile {
        TString::Format("./Outputs/tree_ex_%s_%s_%s_filtered.root", beam.c_str(), target.c_str(), light.c_str())};
    dfFilter.Snapshot("Final_Tree", outfile);
    std::cout << "Saving Final_Tree in " << outfile << '\n';
    // Create and save Ex histos on file
    auto hExSil {dfFilter.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == true; }, {"MergerData"})
                     .Histo1D(HistConfig::Ex, "Ex")};
    hExSil->SetTitle("Ex with silicons");
    auto hExL1 {dfFilter.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"})
                    .Histo1D(HistConfig::Ex, "Ex")};
    hExL1->SetTitle("Ex with L1");
    auto file {std::make_shared<TFile>(outfile.Data(), "update")};
    hExSil->Write("hExSil");
    hExL1->Write("hExL1");
    file->Close();

    // Comparar histogramas
    auto hExBefore = df.Histo1D({"hExBefore", "Excitation Energy before filtering;Ex (MeV);Counts", 100, -5, 10}, "Ex");
    auto hExAfter =
        dfFilter.Histo1D({"hExAfter", "Excitation Energy after filtering;Ex (MeV);Counts", 100, -5, 10}, "Ex");
    if(onlySil)
    {
        hExBefore =
            df.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == true; }, {"MergerData"})
                .Histo1D({"hExBefore", "Excitation Energy before filtering;Ex (MeV);Counts", 100, -5, 10}, "Ex");
        hExAfter = dfFilter.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == true; }, {"MergerData"})
                       .Histo1D({"hExAfter", "Excitation Energy after filtering;Ex (MeV);Counts", 100, -5, 10}, "Ex");
    }

    auto c = new TCanvas("cExFilter", "cExFilter", 800, 600);
    c->Divide(2, 1);
    c->cd(1);
    hExBefore->DrawClone();
    c->cd(2);
    hExAfter->DrawClone();

    auto c2 = new TCanvas("c2ExFilterLog", "c2ExFilterLog", 800, 600);
    c2->Divide(2, 1);
    c2->cd(1);
    hExSil->DrawClone();
    c2->cd(2);
    hExL1->DrawClone();

    // Opcional: guardar eventos rechazados
    // WriteRejectedEvents(infile.Data());
}

#endif
