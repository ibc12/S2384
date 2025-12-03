#ifndef PIPE3_FILTER_H
#define PIPE3_FILTER_H
#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActTPCData.h"
#include "ActSRIM.h"
#include "ActKinematics.h"
#include "ActParticle.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TH1.h"
#include "TMath.h"
#include "TString.h"

#include <fstream>
#include <iostream>
#include <string>

#include "../../PrettyStyle.C"
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
    // PrettyStyle(false);
    bool savePlots {true};
    bool onlySil {true};

    auto infile {TString::Format("./Outputs/tree_ex_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str())};

    ROOT::EnableImplicitMT();
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
                    .Histo1D(HistConfig::ExZoom, "Ex")};
    hExL1->SetTitle("Ex with L1");
    auto file {std::make_shared<TFile>(outfile.Data(), "update")};
    hExSil->Write("hExSil");
    hExL1->Write("hExL1");
    file->Close();
    // Create kin histos for checking
    auto hkinSil {dfFilter.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == true; }, {"MergerData"})
                      .Histo2D(HistConfig::KinPlot, "MergerData.fThetaLight", "EVertex")};
    auto hkinL1 {dfFilter.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"})
                     .Histo2D(HistConfig::KinPlot, "MergerData.fThetaLight", "EVertex")};

    // Comparar histogramas
    auto hExBefore = df.Histo1D({"hExBefore", "Excitation Energy before filtering;E_{x} (MeV);Counts", 100, -5, 10}, "Ex");
    auto hExAfter =
        dfFilter.Histo1D({"hExAfter", "Excitation Energy after filtering;E_{x} (MeV);Counts", 100, -5, 10}, "Ex");
    if(onlySil)
    {
        hExBefore = df.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == true; }, {"MergerData"})
                        .Histo1D({"hExBefore",
                                  TString::Format("Excitation Energy before filtering;E_{x} (MeV);Counts / %.f keV ",
                                                  (10. - (-5.)) / 100 * 1000),
                                  100, -5, 10},
                                 "Ex");
        hExAfter = dfFilter.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == true; }, {"MergerData"})
                       .Histo1D(HistConfig::ExZoom, "Ex");
    }

    // Get theoretical kinematic line for the plots
    ActPhysics::Particle pb {beam};
    ActPhysics::Particle pt {target};
    ActPhysics::Particle pl {light};
    auto* srim {new ActPhysics::SRIM};
    srim->ReadTable(beam, TString::Format("../Calibrations/SRIM/%s_900mb_CF4_95-5.txt", beam.c_str()).Data());
    srim->ReadTable("mylar", TString::Format("../Calibrations/SRIM/%s_Mylar.txt", beam.c_str()).Data());
    // Initial energy
    double initialEnergy {7.558}; // meassured by operators; resolution of 0,19%
    initialEnergy = srim->Slow("mylar", initialEnergy * pb.GetAMU(), 0.0168);
    initialEnergy = srim->Slow(beam, initialEnergy, 60); // 60 mm of gas before the pad plane
    initialEnergy = initialEnergy / pb.GetAMU();         // back to amu units
    ActPhysics::Kinematics kin {pb, pt, pl, initialEnergy * pb.GetAMU()};
    ActPhysics::Kinematics kin1st {pb, pt, pl, initialEnergy * pb.GetAMU(), 1};
    ActPhysics::Kinematics kin2nd {pb, pt, pl, initialEnergy * pb.GetAMU(), 2.2};

    auto c = new TCanvas("cExFilter", "cExFilter", 800, 600);
    c->Divide(2, 1);
    c->cd(1);
    hExBefore->DrawClone();
    c->cd(2);
    hExAfter->DrawClone();

    auto c2 = new TCanvas("c2ExFilter", "c2ExFilter", 800, 600);
    c2->Divide(2, 1);
    c2->cd(1);
    hExSil->DrawClone();
    c2->cd(2);
    hExL1->DrawClone();

    auto c3 = new TCanvas("c3ExFilter", "c3ExFilter", 800, 600);
    c3->Divide(2, 1);
    c3->cd(1);
    hkinSil->DrawClone("colz");
    c3->cd(2);
    hkinL1->DrawClone("colz");

    // Opcional: guardar eventos rechazados
    // WriteRejectedEvents(infile.Data());

    // Save canvases
    if(savePlots)
    {
        // Clone to modify the copy
        auto* hcloneEx = (TH2D*)hExAfter->Clone("hcloneEx");
        auto* hcloneExL1 = (TH2D*)hExL1->Clone("hcloneExL1");
        auto* hcloneKinSil = (TH2D*)hkinSil->Clone("hcloneKinSil");
        auto* hcloneKinL1 = (TH2D*)hkinL1->Clone("hcloneKinL1");

        // Crear canvas temporal solo para guardar
        auto* ctmpEx = new TCanvas("ctmpEx", "", 1600, 1200);
        ctmpEx->cd();
        hcloneEx->DrawClone("hist");
        TString outputEx =
            TString::Format("../Figures/ex_latSil_%s_%s_%s.png", beam.c_str(), target.c_str(), light.c_str());
        ctmpEx->SaveAs(outputEx, "PNG");

        auto* ctmpExL1 = new TCanvas("ctmpExL1", "", 1600, 1200);
        ctmpExL1->cd();
        hcloneExL1->DrawClone("hist");
        TString outputExL1 =
            TString::Format("../Figures/ex_L1_%s_%s_%s.png", beam.c_str(), target.c_str(), light.c_str());
        ctmpExL1->SaveAs(outputExL1, "PNG");

        auto* ctmpKinSil = new TCanvas("ctmpKinSil", "", 1600, 1200);
        ctmpKinSil->cd();
        hcloneKinSil->DrawClone("colz");
        auto* theo {kin.GetKinematicLine3()};
        theo->Draw("same");
        // auto* theo1st {kin1st.GetKinematicLine3()};
        // theo1st->Draw("same");
        // auto* theo2nd {kin2nd.GetKinematicLine3()};
        // theo2nd->Draw("same");
        TString outputKinSil =
            TString::Format("../Figures/kin_latSil_%s_%s_%s.png", beam.c_str(), target.c_str(), light.c_str());
        ctmpKinSil->SaveAs(outputKinSil, "PNG");

        auto* ctmpKinL1 = new TCanvas("ctmpKinL1", "", 1600, 1200);
        ctmpKinL1->cd();
        hcloneKinL1->DrawClone("colz");
        theo->Draw("same");
        // theo1st->Draw("same");
        // theo2nd->Draw("same");
        TString outputKinL1 =
            TString::Format("../Figures/kin_L1_%s_%s_%s.png", beam.c_str(), target.c_str(), light.c_str());
        ctmpKinL1->SaveAs(outputKinL1, "PNG");

        // Guardar
    }
}

#endif
