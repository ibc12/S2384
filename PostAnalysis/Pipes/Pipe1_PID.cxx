#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TH2.h"
#include "TH2D.h"
#include "TString.h"

#include <fstream>
#include <map>
#include <string>

#include "../../PrettyStyle.C"

void Pipe1_PID(const std::string& beam, const std::string& target, const std::string& light)
{
    // PrettyStyle(false);
    bool savePlots = false;
    std::string dataconf {};
    if(beam == "11Li")
        dataconf = "./../configs/data_11Li.conf";
    else if(beam == "7Li")
        dataconf = "./../configs/data_7Li.conf";
    else
        throw std::runtime_error("Beam cannot differ from 11Li or 7Li");

    // Read data
    ActRoot::DataManager dataman {dataconf, ActRoot::ModeType::EMerge};
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    auto chain3 {dataman.GetChain(ActRoot::ModeType::EFilter)};
    auto chain4 {dataman.GetChain(ActRoot::ModeType::EReadTPC)};
    chain->AddFriend(chain2.get());
    chain->AddFriend(chain3.get(), "TPCData");
    chain->AddFriend(chain4.get(), "GETTree");

    // RDataFrame
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame dforigin {*chain};

    // Filter silicon pads
    auto df = dforigin.Filter(
        [](ActRoot::MergerData& m)
        {
            // Mask L0_9
            if(m.fLight.fNs.size())
                if((m.fLight.fLayers.front() == "l0") && (m.fLight.fNs.front() == 9))
                    return false;
            // Mask F0_2
            if(m.fLight.fNs.size())
                if((m.fLight.fLayers.front() == "f0") && (m.fLight.fNs.front() == 2))
                    return false;
            // Mask R0_3
            if(m.fLight.fNs.size())
                if((m.fLight.fLayers.front() == "r0") && (m.fLight.fNs.front() == 3))
                    return false;
            // Mask R0 sils depending on run number
            if(m.fRun > 29 && m.fRun < 35)
            {
                if(m.fRun == 30 || m.fRun == 31 || m.fRun == 33)
                {
                    if(!m.fLight.fLayers.empty() && m.fLight.fLayers.front() == "r0")
                    {
                        if(!m.fLight.fNs.empty() && (m.fLight.fNs.front() == 3 || m.fLight.fNs.front() == 5))
                        {
                            return false;
                        }
                        else
                            return true;
                    }
                    else
                        return true;
                }
                if(m.fRun == 32)
                {
                    if(!m.fLight.fLayers.empty() && m.fLight.fLayers.front() == "r0")
                    {
                        return false;
                    }
                    else
                        return true;
                }
                if(m.fRun == 34)
                {
                    if(!m.fLight.fLayers.empty() &&
                       (m.fLight.fLayers.front() == "r0" || m.fLight.fLayers.front() == "l0"))
                    {
                        return false;
                    }
                    else
                        return true;
                }
                else
                    return true;
            }
            if(m.fRun > 35 && m.fRun < 45)
            {
                if(!m.fLight.fLayers.empty() && m.fLight.fLayers.front() == "f0")
                {
                    if(!m.fLight.fNs.empty() && m.fLight.fNs.front() == 5)
                    {
                        return false;
                    }
                    else
                        return true;
                }
                else
                    return true;
            }
            if(m.fRun == 116 || m.fRun == 117)
            {
                if(!m.fLight.fLayers.empty() && m.fLight.fLayers.front() == "r0")
                {
                    if(!m.fLight.fNs.empty() && m.fLight.fNs.front() == 2)
                    {
                        return false;
                    }
                    else
                        return true;
                }
                else
                    return true;
            }
            else
                return true;
        },
        {"MergerData"});


    // LIGHT particle
    // Define lambda functions
    // 1-> Stopped in first silicon layer
    auto lambdaOne {[](ActRoot::MergerData& m) { return m.fLight.GetNLayers() == 1; }};
    // 2-> In two layers
    auto lambdaTwo {[](ActRoot::MergerData& m)
                    {
                        if(m.fLight.GetNLayers() == 2)
                            return (m.fLight.GetLayer(0) == "f0" && m.fLight.GetLayer(1) == "f1");
                        else
                            return false;
                    }};
    // 3-> Heavy particle
    auto lambdaHeavy {[](ActRoot::MergerData& m) { return m.fHeavy.GetNLayers() == 2; }};
    // 4-> L1
    auto lambdaIsL1 {[](ActRoot::MergerData& mer, ActRoot::ModularData& mod)
                     { return (mod.Get("GATCONF") == 8) && (mer.fLightIdx != -1); }};

    // Fill histograms
    std::map<std::string, ROOT::TThreadedObject<TH2D>> hsgas, hstwo, hszero;
    // Histogram models
    auto hGasSil {new TH2D {"hGasSil", ";E_{Sil} [MeV];#Delta E_{gas} [arb. units]", 450, 0, 70, 600, 0, 3000}};
    auto hTwoSils {new TH2D {"hTwoSils", ";#DeltaE_{0} [MeV];#DeltaE_{1} [MeV]", 500, 0, 80, 400, 0, 35}};
    for(const auto& layer : {"f0", "l0", "r0"})
    {
        hsgas.emplace(layer, *hGasSil);
        hsgas[layer]->SetTitle(TString::Format("%s", layer));
    }
    hstwo.emplace("f0-f1", *hTwoSils);
    hstwo["f0-f1"]->SetTitle("f0-f1");
    for(int s = 0; s < 4; s++)
    {
        hszero.emplace(std::to_string(s), *hTwoSils);
        hszero[std::to_string(s)]->SetTitle(TString::Format("f2_%d vs f3;E_{f3} [MeV];#DeltaE_{f2} [MeV]", s));
    }
    ROOT::TThreadedObject<TH2D> hl1 {"hl1", "L1 PID;Raw TL [au];Q_{total} [au]", 200, 0, 120, 2000, 0, 3e5};
    ROOT::TThreadedObject<TH2D> hl1Gated {"hl1", "L1 PID > 100#circ;Raw TL [au];Q_{total} [au]", 200, 0, 120, 2000, 0,
                                          3e5};
    ROOT::TThreadedObject<TH2D> hl1theta {
        "hl1theta", "L1 #theta;#theta_{L1} [#circ];Q_{total} [au]", 240, 0, 180, 2000, 0, 3e5};
    ROOT::TThreadedObject<TH2D> hl1thetaCorr {
        "hl1thetaCorr", "L1 #thetas;#theta_{Light} [#circ];#theta_{Heavy} [#circ]", 240, 0, 180, 200, 0, 100};

    // Fill them
    df.Foreach(
        [&](ActRoot::MergerData& m, ActRoot::ModularData& mod)
        {
            // L1
            if(lambdaIsL1(m, mod))
            {
                hl1->Fill(m.fLight.fRawTL, m.fLight.fQtotal);
                hl1theta->Fill(m.fThetaLight, m.fLight.fQtotal);
                hl1thetaCorr->Fill(m.fThetaLight, m.fThetaHeavy);
                if(m.fThetaLight > 100)
                    hl1Gated->Fill(m.fLight.fRawTL, m.fLight.fQtotal);
                return;
            }
            // Light
            if(lambdaOne(m)) // Gas-E0 PID
            {
                auto layer {m.fLight.GetLayer(0)};
                if(hsgas.count(layer))
                    hsgas[layer]->Fill(m.fLight.fEs.front(), m.fLight.fQave);
            }
            else if(lambdaTwo(m)) // E0-E1 PID
            {
                hstwo["f0-f1"]->Fill(m.fLight.fEs[0], m.fLight.fEs[1]);
            }
            // Heavy
            if(lambdaHeavy(m))
            {
                auto n {std::to_string(m.fHeavy.fNs[0])};
                if(hszero.count(n))
                    hszero[n]->Fill(m.fHeavy.fEs[1], m.fHeavy.fEs[0]);
            }
        },
        {"MergerData", "ModularData"});

    // If cuts are present, apply them
    ActRoot::CutsManager<std::string> cuts;
    // Gas PID
    cuts.ReadCut("l0", TString::Format("./Cuts/pid_%s_l0_%s.root", light.c_str(), beam.c_str()).Data());
    cuts.ReadCut("r0", TString::Format("./Cuts/pid_%s_r0_%s.root", light.c_str(), beam.c_str()).Data());
    cuts.ReadCut("f0", TString::Format("./Cuts/pid_%s_f0_%s.root", light.c_str(), beam.c_str()).Data());
    cuts.ReadCut("l1", TString::Format("./Cuts/pid_%s_l1_%s.root", light.c_str(), beam.c_str()).Data());

    // Read indivitual cuts for heavy particle --- NOW IN PIPE3
    // if(beam == "11Li" && light == "p")
    // {
    //     for(const auto& heavy : {"11Li", "9Li"}) // these two particles are bound to (d,p) reaction
    //     {
    //         for(int s = 0; s < 4; s++) // one cut per f2 quad pad
    //         {
    //             cuts.ReadCut(TString::Format("%s_f2_%d_f3", heavy, s).Data(),
    //                          TString::Format("./Cuts/pid_%s_f2_%d_%s.root", heavy, s, beam.c_str()).Data());
    //         }
    //     }
    // }

    // Two sils PID
    // cuts.ReadCut("f0-f1", TString::Format("./Cuts/pid_%s_f0_f1_%s.root", light.c_str(), beam.c_str()).Data());
    // Get list of cuts
    auto listOfCuts {cuts.GetListOfKeys()};
    if(listOfCuts.size())
    {
        // Apply PID and save in file
        auto gated {df.Filter(
            [&](ActRoot::MergerData& m, ActRoot::ModularData& mod)
            {
                // L1
                if(lambdaIsL1(m, mod))
                {
                    if(cuts.GetCut("l1"))
                        return cuts.IsInside("l1", m.fLight.fRawTL, m.fLight.fQtotal);
                    else
                        return false;
                }
                // One silicon
                else if(lambdaOne(m))
                {
                    auto layer {m.fLight.GetLayer(0)};
                    if(cuts.GetCut(layer))
                    {
                        // LIGHT particle
                        auto l {cuts.IsInside(layer, m.fLight.fEs[0], m.fLight.fQave)};
                        // HEAVY for (d,p) NOT USED, now in Pipe3
                        // bool h {true};
                        // if(lambdaHeavy(m) && light == "p")
                        // {
                        //     auto n {m.fHeavy.fNs[0]};
                        //     std::set<bool> isInHeavyCuts;
                        //     for(const auto& particle : {"11Li", "9Li"})
                        //     {
                        //         std::string key {TString::Format("%s_f2_%d_f3", particle, n).Data()};
                        //         if(std::find(listOfCuts.begin(), listOfCuts.end(), key) != listOfCuts.end())
                        //             isInHeavyCuts.insert(cuts.IsInside(key, m.fHeavy.fEs[1], m.fHeavy.fEs[0]));
                        //     }
                        //     if(isInHeavyCuts.size())
                        //     {
                        //         if(isInHeavyCuts.find(true) != isInHeavyCuts.end())
                        //         {
                        //             h = true;
                        //         }
                        //         else
                        //             h = false;
                        //     }
                        // }
                        return l;
                    }
                    else
                        return false;
                }
                else if(cuts.GetCut("f0-f1") && lambdaTwo(m)) // PID in fo-f1
                    return cuts.IsInside("f0-f1", m.fLight.fEs[0], m.fLight.fEs[1]);
                else
                    return false;
            },
            {"MergerData", "ModularData"})};
        auto name {TString::Format("./Outputs/tree_pid_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str())};
        std::cout << "Saving PID_Tree in file : " << name << '\n';
        gated.Snapshot("PID_Tree", name.Data());
    }

    // Draw
    auto* c0 {new TCanvas {"c10", "Pipe 1 PID canvas 0"}};
    c0->DivideSquare(6);
    int p {1};
    c0->cd(1);
    for(auto& [layer, h] : hsgas)
    {
        c0->cd(p);
        h.Merge()->DrawClone("colz");
        cuts.DrawCut(layer);
        p++;
    }
    for(auto& [layer, h] : hstwo)
    {
        c0->cd(p);
        h.Merge()->DrawClone("colz");
        p++;
    }
    c0->cd(5);
    hl1.Merge()->DrawClone("colz");
    cuts.DrawCut("l1");
    c0->cd(6);
    hl1theta.Merge()->DrawClone("colz");

    auto* c1 {new TCanvas {"c11", "Pipe1 PID canvas 1"}};
    c1->DivideSquare(hszero.size());
    p = 1;
    for(auto& [s, h] : hszero)
    {
        c1->cd(p);
        h.Merge()->DrawClone("colz");
        for(const auto& particle : {"11Li", "9Li"})
        {
            auto key {TString::Format("%s_f2_%s_f3", particle, s.c_str())};
            cuts.DrawCut(key.Data());
        }
        p++;
    }
    auto c1All {new TCanvas {"c1All", "Pipe1 PID canvas 1 all"}};
    c1All->cd();

    // Obtener el merge del slot 0 como unique_ptr
    auto h0 = hszero["0"].Merge();
    auto h1 = hszero["1"].Merge();
    auto h2 = hszero["2"].Merge();
    auto h3 = hszero["3"].Merge();

    // Crear histograma combinado como copia del slot 0
    auto* hsum = (TH2D*)h0->Clone("hsum");
    hsum->Reset(); // vaciarlo antes de sumar

    // Sumar los histogramas 0 y 2
    hsum->Add(h0.get());
    hsum->Add(h1.get());
    hsum->Add(h2.get());
    hsum->Add(h3.get());

    // Dibujar el histograma combinado
    hsum->Draw("colz");
    // auto c1All {new TCanvas {"c1All", "Pipe1 PID canvas 1 all"}};
    // c1All->cd();
    // p = 0;
    // for(auto& [s, h] : hszero)
    // {
    //     if(p == 0)
    //         h.Merge()->DrawClone("colz");
    //     else if(p == 3)
    //         continue;
    //     else if(p == 1)
    //         continue;
    //     else if(p == 2)
    //         h.Merge()->DrawClone("SAME");
    //     for(const auto& particle : {"11Li", "9Li"})
    //     {
    //         auto key {TString::Format("%s_f2_%s_f3", particle, s.c_str())};
    //         cuts.DrawCut(key.Data());
    //     }
    //     p++;
    // }

    auto* c2 {new TCanvas {"c12", "Pipe1 PID canvas 2"}};
    c2->DivideSquare(4);
    c2->cd(1);
    hl1.GetAtSlot(0)->DrawClone("colz");
    c2->cd(2);
    hl1theta.GetAtSlot(0)->DrawClone("colz");
    c2->cd(3);
    hl1thetaCorr.Merge()->DrawClone("colz");
    auto* gtheo {ActPhysics::Kinematics(TString::Format("%s(d,p)@82.5", beam.c_str()).Data()).GetTheta3vs4Line()};
    gtheo->SetLineColor(46);
    gtheo->Draw("l");
    c2->cd(4);
    hl1Gated.Merge()->DrawClone("colz");

    // Save important canvases
    if(savePlots)
    {
        // Get the histograms
        auto hptrL0 = hsgas.at("l0").Merge();
        auto hptrR0 = hsgas.at("r0").Merge();
        auto hptrF2F3 = hszero.at("0").Merge();

        // Clone to modify axis ranges
        auto* hmodL0 = (TH2D*)hptrL0->Clone("hmodL0");
        auto* hmodR0 = (TH2D*)hptrR0->Clone("hmodR0");
        auto* hmodF2F3 = (TH2D*)hptrF2F3->Clone("hmodF2F3");

        hmodL0->GetXaxis()->SetRangeUser(0, 16);
        hmodL0->GetYaxis()->SetRangeUser(0, 700);

        hmodR0->GetXaxis()->SetRangeUser(0, 16);
        hmodR0->GetYaxis()->SetRangeUser(0, 700);

        hmodF2F3->GetXaxis()->SetRangeUser(0, 80);
        hmodF2F3->GetYaxis()->SetRangeUser(0, 20);

        auto* ctmpL0 = new TCanvas("ctmpL0", "", 1600, 1200);
        ctmpL0->cd();
        hmodL0->Draw("colz");
        ctmpL0->Update();
        TString outputL0 = TString::Format("../Figures/pid_l0_%s.png", beam.c_str());
        ctmpL0->SaveAs(outputL0);

        auto* ctmpR0 = new TCanvas("ctmpR0", "", 1600, 1200);
        ctmpR0->cd();
        hmodR0->Draw("colz");
        ctmpR0->Update();
        TString outputR0 = TString::Format("../Figures/pid_r0_%s.png", beam.c_str());
        ;
        ctmpR0->SaveAs(outputR0);

        auto* ctmpF2F3 = new TCanvas("ctmpF2F3", "", 1600, 1200);
        ctmpF2F3->cd();
        hmodF2F3->Draw("colz");
        ctmpF2F3->Update();
        TString outputF2F3 = TString::Format("../Figures/pid_f2f3_%s.png", beam.c_str());
        ctmpF2F3->SaveAs(outputF2F3);
    }

    // ======================================
    // DEBUG: f0 / l0 / r0 con HEAVY coincidente
    // ======================================

    // DataFrame: solo eventos
    //  - 1 silicio (gas–E0 PID)
    //  - layer = f0, l0 o r0
    //  - con detección de heavy (lambdaHeavy ya existente)
    // auto df_f0l0r0_heavy = df.Filter(
    //    [&](ActRoot::MergerData& m)
    //    {
    //        // Heavy detectado
    //        if(!lambdaHeavy(m))
    //            return false;
    //
    //        // Un solo silicio
    //        if(m.fLight.GetNLayers() != 1)
    //            return false;
    //
    //        // Solo f0, l0 o r0
    //        auto layer = m.fLight.GetLayer(0);
    //        return (layer == "f0" || layer == "l0" || layer == "r0");
    //    },
    //    {"MergerData"});
    //
    //// Histogramas debug (gas–sil) SOLO con heavy
    // std::map<std::string, ROOT::TThreadedObject<TH2D>> hsgasHeavy;
    //
    // for(const auto& layer : {"f0", "l0", "r0"})
    //{
    //    hsgasHeavy.emplace(layer, TH2D(TString::Format("hGasSil_%s_heavy", layer),
    //                                   TString::Format("%s + heavy;E_{Sil} [MeV];#Delta E_{gas} [au]", layer), 450, 0,
    //                                   70, 600, 0, 3000));
    //}
    //
    //// Llenado
    // cuts.ReadCut("alfas", TString::Format("./Cuts/pid_all_alfas_7Li.root").Data());
    // df_f0l0r0_heavy.Foreach(
    //     [&](ActRoot::MergerData& m)
    //     {
    //         if(cuts.IsInside("alfas", m.fHeavy.fEs[1], m.fHeavy.fEs[0]))
    //         {
    //             auto layer = m.fLight.GetLayer(0);
    //             hsgasHeavy[layer]->Fill(m.fLight.fEs.front(), m.fLight.fQave);
    //         }
    //     },
    //     {"MergerData"});
    //
    //// Canvas debug: superposición
    // auto* cDebug = new TCanvas("cDebug_f0l0r0_heavy", "DEBUG f0/l0/r0 with heavy", 1800, 600);
    // cDebug->Divide(3, 1);
    //
    //// f0
    // cDebug->cd(1);
    // auto h_f0_all = hsgas["f0"].Merge();
    // auto h_f0_heavy = hsgasHeavy["f0"].Merge();
    //
    // h_f0_all->SetTitle("f0 PID (all) + heavy overlay");
    // h_f0_all->DrawClone("colz");
    //
    // h_f0_heavy->SetLineColor(kRed);
    // h_f0_heavy->SetLineWidth(3);
    // h_f0_heavy->DrawClone("cont3 SAME");
    //
    // cuts.DrawCut("f0");
    //
    //// l0
    // cDebug->cd(2);
    // auto h_l0_all = hsgas["l0"].Merge();
    // auto h_l0_heavy = hsgasHeavy["l0"].Merge();
    //
    // h_l0_all->SetTitle("l0 PID (all) + heavy overlay");
    // h_l0_all->DrawClone("colz");
    //
    // h_l0_heavy->SetLineColor(kRed);
    // h_l0_heavy->SetLineWidth(3);
    // h_l0_heavy->DrawClone("cont3 SAME");
    //
    // cuts.DrawCut("l0");
    //
    //// r0
    // cDebug->cd(3);
    // auto h_r0_all = hsgas["r0"].Merge();
    // auto h_r0_heavy = hsgasHeavy["r0"].Merge();
    //
    // h_r0_all->SetTitle("r0 PID (all) + heavy overlay");
    // h_r0_all->DrawClone("colz");
    //
    // h_r0_heavy->SetLineColor(kRed);
    // h_r0_heavy->SetLineWidth(3);
    // h_r0_heavy->DrawClone("cont3 SAME");
    //
    // cuts.DrawCut("r0");
    //
    //// Save some events if needed
    // cuts.ReadCut("f0_z1", TString::Format("./Cuts/cut_z1Particles_11Li.root").Data());
    // auto dfOut = df.Filter(
    //     [&](ActRoot::MergerData& m)
    //     {
    //         if(m.fLight.GetLayer(0) == "f0" && m.fLight.GetNLayers() == 1)
    //             return (cuts.IsInside("f0_z1", m.fLight.fEs[0], m.fLight.fQave));
    //         else
    //             return false;
    //     },
    //     {"MergerData"});
    // std::ofstream out("./Outputs/pid_events_z1Front_11Li.dat");
    // dfOut.Foreach([&](ActRoot::MergerData& m) { m.Stream(out); }, {"MergerData"});
    // out.close();
    // cuts.ReadCut("f0f1", TString::Format("./cut_f0f1_events.root").Data());
    // auto dfOut = df.Filter(
    //     [&](ActRoot::MergerData& m)
    //     {
    //         if(m.fLight.GetLayer(0) == "f0" && m.fLight.GetLayer(1) == "f1" && m.fLight.GetNLayers() == 2)
    //             return (cuts.IsInside("f0f1", m.fLight.fEs[0], m.fLight.fEs[1]));
    //         else
    //             return false;
    //     },
    //     {"MergerData"});
    //std::ofstream out("./Outputs/pid_events_f0f1_11Li.dat");
    //dfOut.Foreach([&](ActRoot::MergerData& m) { m.Stream(out); }, {"MergerData"});
    //out.close();
}
