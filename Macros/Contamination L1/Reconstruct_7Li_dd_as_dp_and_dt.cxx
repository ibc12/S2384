#include "ActKinematicGenerator.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActParticle.h"
#include "ActSRIM.h"

#include "ROOT/RDF/HistoModels.hxx"
#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TGraph.h"
#include "TH2.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPad.h"
#include "TROOT.h"
#include "TString.h"

#include "FitInterface.h"
#include "FitModel.h"
#include "FitUtils.h"
#include "Interpolators.h"

#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "../../Fits/Histos.h"
#include "../../PostAnalysis/HistConfig.h"
#include "../../PrettyStyle.C"

void Reconstruct_7Li_dd_as_dp_and_dt()
{
    PrettyStyle(false, true);
    ROOT::EnableImplicitMT();

    // ------------------------------------------------------------------
    // 1. Analysis: read tree and filters
    //    (los arboles ya salen de Pipe2 con reconstruccion dd:
    //     EVertex, EBeam, ThetaCM, Ex ya existen)
    // ------------------------------------------------------------------
    ROOT::RDataFrame df_ini {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_7Li_d_d_filtered.root"};
    auto df {df_ini.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"})};
    ROOT::RDataFrame df_extended_ini {"Final_Tree", "./Inputs/tree_ex_F_7Li_d_d_filtered_extendedUpCutInQtot.root"};
    auto df_extended {
        df_extended_ini.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"})};

    // ------------------------------------------------------------------
    // 2. Prepare all the kinematics, dd, dp, dt
    //    (la energia del haz se sobreescribe evento a evento con EBeam)
    // ------------------------------------------------------------------
    const std::string beam {"7Li"};
    ActPhysics::Particle pb {beam};
    ActPhysics::Particle pt {"d"};
    ActPhysics::Particle pl_d {"d"};
    ActPhysics::Particle pl_t {"t"};
    ActPhysics::Particle pl_p {"p"};
    double initialEnergy {7.225};
    ActPhysics::Kinematics kin_dd {pb, pt, pl_d, initialEnergy * pb.GetAMU()};
    ActPhysics::Kinematics kin_dp {pb, pt, pl_p, initialEnergy * pb.GetAMU()};
    ActPhysics::Kinematics kin_dt {pb, pt, pl_t, initialEnergy * pb.GetAMU()};

    // ------------------------------------------------------------------
    // 3. Declare all srim objects to reconstruct EVertex
    //    (ya no hace falta la tabla del haz ni la de mylar: EBeam ya esta calculado)
    // ------------------------------------------------------------------
    auto* srim {new ActPhysics::SRIM};
    srim->ReadTable("p", "../../Calibrations/SRIM/1H_900mb_CF4_95-5.txt");
    srim->ReadTable("t", "../../Calibrations/SRIM/3H_900mb_CF4_95-5.txt");

    // ------------------------------------------------------------------
    // 4. Defines de EVertex, ThetaCM y Ex para las hipotesis p y t
    //    (d = alias de las columnas ya existentes)
    // ------------------------------------------------------------------
    struct Hypo
    {
        std::string suffix; // p, t
        std::string srimKey;
        ActPhysics::Kinematics kin;
    };
    std::vector<Hypo> hypos {{"p", "p", kin_dp}, {"t", "t", kin_dt}};

    auto addVariables = [&](ROOT::RDF::RNode node) -> ROOT::RDF::RNode
    {
        const auto nslots {node.GetNSlots()};

        // Hipotesis dd: ya calculada en Pipe2 -> solo alias
        node = node.Alias("EVertex_d", "EVertex");
        node = node.Alias("ThetaCM_d", "ThetaCM");
        node = node.Alias("Ex_d", "Ex");

        for(auto& h : hypos)
        {
            const std::string key {h.srimKey};
            const std::string sfx {h.suffix};

            // EVertex_x: misma logica que Pipe2 pero con la tabla SRIM de la particula supuesta
            node = node.Define("EVertex_" + sfx,
                               [srim, key](const ActRoot::MergerData& d)
                               {
                                   double ret {};
                                   if(d.fLight.IsFilled())
                                       ret = srim->EvalInitialEnergy(key, d.fLight.fEs.front(), d.fLight.fTL);
                                   else // L1 trigger
                                       ret = srim->EvalEnergy(key, d.fLight.fTL);
                                   return ret;
                               },
                               {"MergerData"});

            // Un objeto Kinematics por slot (se cambia la energia del haz en cada entrada)
            auto vkins {std::make_shared<std::vector<ActPhysics::Kinematics>>(nslots, h.kin)};

            // ThetaCM_x (usa la columna EBeam ya existente)
            node = node.DefineSlot(
                "ThetaCM_" + sfx,
                [vkins](unsigned int slot, const ActRoot::MergerData& d, double EVertex, double EBeam)
                {
                    auto& k {(*vkins)[slot]};
                    k.SetBeamEnergy(EBeam);
                    return k.ReconstructTheta3CMFromLab(EVertex, d.fThetaLight * TMath::DegToRad()) * TMath::RadToDeg();
                },
                {"MergerData", "EVertex_" + sfx, "EBeam"});

            // Ex_x
            node =
                node.DefineSlot("Ex_" + sfx,
                                [vkins](unsigned int slot, const ActRoot::MergerData& d, double EVertex, double EBeam)
                                {
                                    auto& k {(*vkins)[slot]};
                                    k.SetBeamEnergy(EBeam);
                                    return k.ReconstructExcitationEnergy(EVertex, d.fThetaLight * TMath::DegToRad());
                                },
                                {"MergerData", "EVertex_" + sfx, "EBeam"});
        }
        return node;
    };

    auto def {addVariables(df)};
    auto def_ext {addVariables(df_extended)};
    auto def_ext_noElastic {def_ext.Filter([](double ex) { return ex > 1; }, {"Ex"})};

    // ------------------------------------------------------------------
    // 5. Histogramas (solo eventos L1) -> con el arbol EXTENDIDO
    // ------------------------------------------------------------------
    const ROOT::RDF::TH1DModel Ex_modified {
        "hEx", TString::Format("Excitation energy;E_{x} [MeV];Counts / %.f keV", (12. - (-8.)) / 125 * 1e3), 125, -8,
        12};
    const ROOT::RDF::TH2DModel KinEl_modified {"hKinEl", "Kinematics;#theta_{Lab} [#circ];E_{Vertex} [MeV]", 600, 0, 180, 100, 0, 5};
    const ROOT::RDF::TH2DModel KinCM_modified {"hKinCM", "CM kinematics;#theta_{CM} [#circ];E_{Vertex} [MeV]", 400, 0, 60, 100, 0, 5};

    auto hEx_d_all {def_ext.Histo1D(Ex_modified, "Ex_d")};
    auto hEx_d {def_ext_noElastic.Histo1D(Ex_modified, "Ex_d")};
    auto hEx_p {def_ext_noElastic.Histo1D(Ex_modified, "Ex_p")};
    auto hEx_t {def_ext_noElastic.Histo1D(Ex_modified, "Ex_t")};
    hEx_d->SetTitle("Ex as dd (L1, extended)");
    hEx_p->SetTitle("Ex as dp (L1, extended)");
    hEx_t->SetTitle("Ex as dt (L1, extended)");

    auto hKinCM_d {def_ext_noElastic.Histo2D(KinCM_modified, "ThetaCM_d", "EVertex_d")};
    auto hKinCM_p {def_ext_noElastic.Histo2D(KinCM_modified, "ThetaCM_p", "EVertex_p")};
    auto hKinCM_t {def_ext_noElastic.Histo2D(KinCM_modified, "ThetaCM_t", "EVertex_t")};

    auto hKin_d {def_ext_noElastic.Histo2D(KinEl_modified, "fThetaLight", "EVertex_d")};
    auto hKin_p {def_ext_noElastic.Histo2D(KinEl_modified, "fThetaLight", "EVertex_p")};
    auto hKin_t {def_ext_noElastic.Histo2D(KinEl_modified, "fThetaLight", "EVertex_t")};

    // Ex como dd con el arbol normal (solo para la comparacion)
    auto hEx_d_nom {def.Histo1D(Ex_modified, "Ex_d")};

    // ------------------------------------------------------------------
    // 6. Guardar arboles con las nuevas columnas
    // ------------------------------------------------------------------
    // auto outfile {TString::Format("./Outputs/tree_ex_%s_dd_as_dp_dt.root", beam.c_str())};
    // def.Snapshot("Final_Tree", outfile.Data());
    // std::cout << "Saving Final_Tree in " << outfile << '\n';
    //
    // auto outfile_ext {TString::Format("./Outputs/tree_ex_%s_dd_as_dp_dt_extended.root", beam.c_str())};
    // def_ext.Snapshot("Final_Tree", outfile_ext.Data());
    // std::cout << "Saving Final_Tree (extended) in " << outfile_ext << '\n';

    // ------------------------------------------------------------------
    // 7. Dibujar
    // ------------------------------------------------------------------
    auto* c1 {new TCanvas("c1", "dd as dd / dp / dt (extended)")};
    c1->DivideSquare(9);
    c1->cd(1);
    hKin_d->DrawClone("colz");
    kin_dd.GetKinematicLine3()->Draw("same");
    c1->cd(2);
    hKin_p->DrawClone("colz");
    kin_dp.GetKinematicLine3()->Draw("same");
    c1->cd(3);
    hKin_t->DrawClone("colz");
    kin_dt.GetKinematicLine3()->Draw("same");
    c1->cd(4);
    hKinCM_d->DrawClone("colz");
    c1->cd(5);
    hKinCM_p->DrawClone("colz");
    c1->cd(6);
    hKinCM_t->DrawClone("colz");
    c1->cd(7);
    hEx_d->DrawClone();
    c1->cd(8);
    hEx_p->DrawClone();
    c1->cd(9);
    hEx_t->DrawClone();

    // ------------------------------------------------------------------
    // 8. Comparacion Ex (dd): normal vs extendido
    // ------------------------------------------------------------------
    auto* hNom {static_cast<TH1D*>(hEx_d_nom->Clone("hEx_d_nom_cmp"))};
    auto* hExt {static_cast<TH1D*>(hEx_d_all->Clone("hEx_d_ext_cmp"))};
    hNom->SetLineColor(kBlue);
    hExt->SetLineColor(kRed);
    hNom->SetTitle("Ex as dd (L1): normal vs extended");
    hExt->SetTitle("Ex as dd (L1): normal vs extended");

    // Version normalizada a area (para comparar forma independientemente de la estadistica)
    auto* hNomN {static_cast<TH1D*>(hNom->Clone("hEx_d_nom_cmp_norm"))};
    auto* hExtN {static_cast<TH1D*>(hExt->Clone("hEx_d_ext_cmp_norm"))};
    if(hNomN->Integral() > 0)
        hNomN->Scale(1. / hNomN->Integral());
    if(hExtN->Integral() > 0)
        hExtN->Scale(1. / hExtN->Integral());
    hNomN->SetTitle("Ex as dd (L1): area normalized");
    hExtN->SetTitle("Ex as dd (L1): area normalized");

    auto* c2 {new TCanvas("c2", "Ex dd: normal vs extended")};
    c2->DivideSquare(2);

    c2->cd(1);
    hExt->SetMaximum(1.1 * std::max(hNom->GetMaximum(), hExt->GetMaximum()));
    hExt->DrawClone("hist");
    hNom->DrawClone("hist same");
    auto* leg1 {new TLegend(0.6, 0.72, 0.89, 0.89)};
    leg1->AddEntry(hNom, TString::Format("Normal (%.0f)", hNom->GetEntries()), "l");
    leg1->AddEntry(hExt, TString::Format("Extended (%.0f)", hExt->GetEntries()), "l");
    leg1->Draw();

    c2->cd(2);
    hExtN->SetMaximum(1.1 * std::max(hNomN->GetMaximum(), hExtN->GetMaximum()));
    hExtN->DrawClone("hist");
    hNomN->DrawClone("hist same");
    auto* leg2 {new TLegend(0.6, 0.72, 0.89, 0.89)};
    leg2->AddEntry(hNomN, "Normal", "l");
    leg2->AddEntry(hExtN, "Extended", "l");
    leg2->Draw();
}