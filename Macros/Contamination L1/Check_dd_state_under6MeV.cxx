#include "ActKinematicGenerator.h"
#include "ActKinematics.h"
#include "ActMergerData.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TGraph.h"
#include "TH2.h"
#include "TLegend.h"
#include "TPad.h"
#include "TROOT.h"

#include "FitInterface.h"
#include "FitModel.h"
#include "FitUtils.h"
#include "Interpolators.h"

#include <fstream>
#include <iostream>
#include <mutex>
#include <string>
#include <vector>

#include "../../Fits/Histos.h"
#include "../../PostAnalysis/HistConfig.h"
#include "../../PrettyStyle.C"

void Check_dd_state_under6MeV()
{
    PrettyStyle(false, true);
    ROOT::EnableImplicitMT();

    // Ex windows for the selections
    const double exMin6 {5.4};
    const double exMax6 {6.0};
    const double exMin24 {2.0};
    const double exMax24 {4.0};

    // ------------------------------------------------------------------
    // 1. Analysis: read tree and filters
    // ------------------------------------------------------------------
    ROOT::RDataFrame df {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_7Li_d_d_filtered.root"};

    // == false is for L1 events
    auto def {df.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"})};

    // Events in the Ex windows
    auto def_6MeV {def.Filter([=](double ex) { return ex < exMax6 && ex > exMin6; }, {"Ex"})};
    auto def_2_4MeV {def.Filter([=](double ex) { return ex < exMax24 && ex > exMin24; }, {"Ex"})};

    // ------------------------------------------------------------------
    // 2. Book all histograms (event loop not triggered yet)
    // ------------------------------------------------------------------
    // Ex
    auto hEx {def.Histo1D(S2384Fit::Exdd_7Li, "Ex")};

    // PID: Qtotal vs RawTL
    auto hQtot_RawTL {
        def.Histo2D({"hQtot_RawTL", "Qtotal vs RawTL;RawTL [mm];Qtotal [a.u.]", 120, 0, 200, 1000, 0, 3e5},
                    "MergerData.fLight.fRawTL", "MergerData.fLight.fQtotal")};
    auto hQtot_RawTL_6MeV {def_6MeV.Histo2D(
        {"hQtot_RawTL_6MeV", "Qtotal vs RawTL for Ex 5.4-6 MeV;RawTL [mm];Qtotal [a.u.]", 120, 0, 200, 1000, 0, 3e5},
        "MergerData.fLight.fRawTL", "MergerData.fLight.fQtotal")};
    auto hQtot_RawTL_2_4MeV {def_2_4MeV.Histo2D(
        {"hQtot_RawTL_2_4MeV", "Qtotal vs RawTL for Ex 2-4 MeV;RawTL [mm];Qtotal [a.u.]", 120, 0, 200, 1000, 0, 3e5},
        "MergerData.fLight.fRawTL", "MergerData.fLight.fQtotal")};

    // Kinematics: EVertex vs theta
    auto hKin {def.Histo2D({"hKin", "Kinematics;#theta_{Lab} [deg];E_{Lab} [MeV]", 180, 0, 180, 100, 0, 5},
                           "MergerData.fThetaLight", "EVertex")};
    auto hKin_6MeV {def_6MeV.Histo2D(
        {"hKin_6MeV", "Kinematics for Ex 5.4-6 MeV;#theta_{Lab} [deg];E_{Lab} [MeV]", 180, 0, 180, 100, 0, 5},
        "MergerData.fThetaLight", "EVertex")};
    auto hKin_2_4MeV {def_2_4MeV.Histo2D(
        {"hKin_2_4MeV", "Kinematics for Ex 2-4 MeV;#theta_{Lab} [deg];E_{Lab} [MeV]", 180, 0, 180, 100, 0, 5},
        "MergerData.fThetaLight", "EVertex")};

    // Ex vs thetaLab and thetaCM
    auto hEx_Theta {
        def.Histo2D({"hEx_Theta", "Ex vs #theta_{Lab};Ex [MeV];#theta_{Lab} [deg]", 100, -2, 10, 180, 0, 180}, "Ex",
                    "MergerData.fThetaLight")};
    auto hEx_ThetaCM {def.Histo2D(
        {"hEx_ThetaCM", "Ex vs #theta_{CM};Ex [MeV];#theta_{CM} [deg]", 100, -2, 10, 180, 0, 180}, "Ex", "ThetaCM")};

    // Kinematic objects to future plot
    ActPhysics::Particle pb {"7Li"};
    ActPhysics::Particle pt {"d"};
    ActPhysics::Particle pl_d {"d"};
    ActPhysics::Particle pl_t {"t"};
    ActPhysics::Particle pl_p {"p"};
    double initialEnergy {7.225}; // meassured by operators; resolution of 0,19%
    ActPhysics::Kinematics kin_dd_5_8 {pb, pt, pl_d, initialEnergy * pb.GetAMU(), 5.8};
    ActPhysics::Kinematics kin_dt {pb, pt, pl_t, initialEnergy * pb.GetAMU()};

    // ------------------------------------------------------------------
    // 3. Raw pointers (the first GetPtr triggers ONE event loop that fills
    //    all the histograms booked above)
    // ------------------------------------------------------------------
    auto* pEx {hEx.GetPtr()};
    auto* pQtot_RawTL {hQtot_RawTL.GetPtr()};
    auto* pQtot_RawTL_6MeV {hQtot_RawTL_6MeV.GetPtr()};
    auto* pQtot_RawTL_2_4MeV {hQtot_RawTL_2_4MeV.GetPtr()};
    auto* pKin {hKin.GetPtr()};
    auto* pKin_6MeV {hKin_6MeV.GetPtr()};
    auto* pKin_2_4MeV {hKin_2_4MeV.GetPtr()};
    auto* pEx_Theta {hEx_Theta.GetPtr()};
    auto* pEx_ThetaCM {hEx_ThetaCM.GetPtr()};

    pEx->SetAxisRange(1, 9, "X");
    std::cout << "Entries in Ex " << exMin6 << "-" << exMax6 << " MeV: " << pQtot_RawTL_6MeV->GetEntries() << '\n';
    std::cout << "Entries in Ex " << exMin24 << "-" << exMax24 << " MeV: " << pQtot_RawTL_2_4MeV->GetEntries() << '\n';

    // ------------------------------------------------------------------
    // 4. TGraphs with one point per filled bin (to overlay on the total PID)
    // ------------------------------------------------------------------
    auto makeGraph {[](TH2* h, int color)
                    {
                        auto* g {new TGraph {}};
                        for(int ix = 1; ix <= h->GetNbinsX(); ix++)
                            for(int iy = 1; iy <= h->GetNbinsY(); iy++)
                                if(h->GetBinContent(ix, iy) > 0)
                                    g->AddPoint(h->GetXaxis()->GetBinCenter(ix), h->GetYaxis()->GetBinCenter(iy));
                        g->SetMarkerStyle(20);
                        g->SetMarkerSize(0.6);
                        g->SetMarkerColor(color);
                        return g;
                    }};
    auto* g6 {makeGraph(pQtot_RawTL_6MeV, kRed)};
    auto* g24 {makeGraph(pQtot_RawTL_2_4MeV, kBlue)};

    // ------------------------------------------------------------------
    // 5. Draw
    // ------------------------------------------------------------------
    // Ex
    auto* c_ex {new TCanvas {"c_Ex", "Ex canvas"}};
    c_ex->cd();
    pEx->DrawClone();

    // PID: total + selections on top
    auto* c_qtot {new TCanvas {"c_Qtot_RawTL", "Qtotal vs RawTL canvas", 1500, 500}};
    c_qtot->Divide(3, 1);

    // Pad 1: total PID + Ex 5.4-6 MeV
    c_qtot->cd(1);
    gPad->SetLogz();
    pQtot_RawTL->DrawClone("colz");
    g6->DrawClone("P"); // no "A": uses the axes of the histogram already drawn
    auto* leg6 {new TLegend {0.55, 0.78, 0.88, 0.88}};
    leg6->AddEntry(g6, Form("Ex %.1f-%.1f MeV", exMin6, exMax6), "p");
    leg6->Draw();

    // Pad 2: total PID + Ex 2-4 MeV
    c_qtot->cd(2);
    gPad->SetLogz();
    pQtot_RawTL->DrawClone("colz");
    g24->DrawClone("P");
    auto* leg24 {new TLegend {0.55, 0.78, 0.88, 0.88}};
    leg24->AddEntry(g24, Form("Ex %.1f-%.1f MeV", exMin24, exMax24), "p");
    leg24->Draw();

    // Pad 3: total PID + both selections
    c_qtot->cd(3);
    gPad->SetLogz();
    pQtot_RawTL->DrawClone("colz");
    g24->DrawClone("P");
    g6->DrawClone("P");
    auto* legBoth {new TLegend {0.55, 0.74, 0.88, 0.88}};
    legBoth->AddEntry(g6, Form("Ex %.1f-%.1f MeV", exMin6, exMax6), "p");
    legBoth->AddEntry(g24, Form("Ex %.1f-%.1f MeV", exMin24, exMax24), "p");
    legBoth->Draw();

    // PID: each selection alone
    auto* c_qtot_sel {new TCanvas {"c_Qtot_RawTL_sel", "Qtotal vs RawTL selections", 1000, 500}};
    c_qtot_sel->Divide(2, 1);
    c_qtot_sel->cd(1);
    pQtot_RawTL_6MeV->DrawClone("colz");
    c_qtot_sel->cd(2);
    pQtot_RawTL_2_4MeV->DrawClone("colz");

    // Kinematics
    auto* c_kin {new TCanvas {"c_Kin", "Kinematics canvas", 1500, 500}};
    c_kin->Divide(3, 1);
    c_kin->cd(1);
    pKin->DrawClone("colz");
    auto line_dd = kin_dd_5_8.GetKinematicLine3();
    line_dd->SetLineColor(kRed);
    line_dd->Draw("same");
    auto line_dt = kin_dt.GetKinematicLine3();
    line_dt->SetLineColor(kBlue);
    line_dt->Draw("same");
    c_kin->cd(2);
    pKin_6MeV->DrawClone("colz");
    c_kin->cd(3);
    pKin_2_4MeV->DrawClone("colz");

    // Ex vs theta
    auto* c_ex_theta {new TCanvas {"c_Ex_Theta", "Ex vs theta canvas"}};
    c_ex_theta->Divide(2, 1);
    c_ex_theta->cd(1);
    pEx_ThetaCM->DrawClone("colz");
    c_ex_theta->cd(2);
    pEx_Theta->DrawClone("colz");

    // ------------------------------------------------------------------
    // 6. Save events to visualize them by hand in actroot
    //    (mutex because Foreach runs in parallel with ImplicitMT and
    //     writing to a single ofstream from several threads is not safe)
    // ------------------------------------------------------------------
    // std::mutex mtx;
    // std::ofstream out("./Outputs/state_between2_4MeV_7Li_dd.dat");
    // def_2_4MeV.Foreach(
    //     [&](ActRoot::MergerData& m)
    //     {
    //         std::lock_guard<std::mutex> lock {mtx};
    //         m.Stream(out);
    //     },
    //     {"MergerData"});
    // out.close();

    // std::ofstream out6("./Outputs/state_under6MeV_7Li_dd.dat");
    // def_6MeV.Foreach(
    //     [&](ActRoot::MergerData& m)
    //     {
    //         std::lock_guard<std::mutex> lock {mtx};
    //         m.Stream(out6);
    //     },
    //     {"MergerData"});
    // out6.close();
}