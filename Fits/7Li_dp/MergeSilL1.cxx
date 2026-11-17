#include "ActKinematics.h"
#include "ActMergerData.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TColor.h"
#include "TGraphErrors.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"

#include "AngComparator.h"
#include "AngDifferentialXS.h"
#include "AngFitter.h"
#include "AngGlobals.h"
#include "AngIntervals.h"
#include "FitInterface.h"
#include "Interpolators.h"
#include "PhysExperiment.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "../../PrettyStyle.C"
#include "../Histos.h"

// Lee un fichero con columnas x y sy y devuelve un TGraphErrors (vacio si no existe)
TGraphErrors* ReadXS(const std::string& path)
{
    auto* g {new TGraphErrors()};
    std::ifstream in(path);
    if(!in.is_open())
    {
        std::cout << "[MergeSilL1] WARNING: no se puede abrir " << path << '\n';
        return g;
    }
    double x, y, sy;
    while(in >> x >> y >> sy)
    {
        int n = g->GetN();
        g->SetPoint(n, x, y);
        g->SetPointError(n, 0., sy); // ex = 0, ey = sy
    }
    return g;
}

// Junta dos graficos: primero L1 y luego sil (para mantener el orden)
TGraphErrors* MergeGraphs(TGraphErrors* gL1, TGraphErrors* gSil)
{
    auto* g {new TGraphErrors()};
    for(auto* src : {gL1, gSil})
    {
        for(int i = 0; i < src->GetN(); ++i)
        {
            double x, y;
            src->GetPoint(i, x, y);
            int n = g->GetN();
            g->SetPoint(n, x, y);
            g->SetPointError(n, src->GetErrorX(i), src->GetErrorY(i));
        }
    }
    return g;
}

// Pipeline para los picos v0..v8: lectura, fusion, comparacion con p-wave y d-wave
void ProcessUnbound(const std::string& peak)
{
    const std::string silPath {"./Outputs/xs/" + peak + "_xs_lDependent.dat"};
    const std::string l1Path {"../7Li_dp_testL1/Outputs/xs/" + peak + "_xs_lDependent.dat"};
    const std::string pWave {"./Inputs/Preliminary_xs_p_d_waves_unbound/p_wave/fresco/fort.202"};
    const std::string dWave {"./Inputs/Preliminary_xs_p_d_waves_unbound/d_wave/fresco/fort.202"};

    auto* gSil {ReadXS(silPath)};
    auto* gL1 {ReadXS(l1Path)};
    if(gSil->GetN() == 0 && gL1->GetN() == 0)
    {
        std::cout << "[MergeSilL1] Sin datos para " << peak << ", se salta\n";
        return;
    }
    auto* gMerged {MergeGraphs(gL1, gSil)};

    // Se crean en el heap para que no se destruyan al salir de la funcion
    // (igual que quedaban vivos al final de MergeSilL1 en la version original)
    auto* comp {new Angular::Comparator("Merged " + peak, gMerged)};
    comp->Add("p-wave", pWave);
    comp->Add("d-wave", dWave);
    comp->Fit();
    comp->Draw("", true);

    auto* compSil {new Angular::Comparator("Sil " + peak, gSil)};
    compSil->Add("p-wave", pWave);
    compSil->Add("d-wave", dWave);
    compSil->Fit();
    compSil->Draw("", true);

    auto* compL1 {new Angular::Comparator("L1 " + peak, gL1)};
    compL1->Add("p-wave", pWave);
    compL1->Add("d-wave", dWave);
    compL1->Fit();
    compL1->Draw("", true);
}

void MergeSilL1()
{
    PrettyStyle();

    // ---------- g0 (ground state) ----------
    auto* g_sil_xs {ReadXS("./Outputs/xs/g0_xs_lDependent.dat")};
    auto* g_L1_xs {ReadXS("../7Li_dp_testL1/Outputs/xs/g0_xs_lDependent.dat")};
    auto* g_merged {MergeGraphs(g_L1_xs, g_sil_xs)};

    // ---------- g1 (1st excited state) ----------
    auto* g_sil_xs_1st {ReadXS("./Outputs/xs/g1_xs_lDependent.dat")};
    auto* g_L1_xs_1st {ReadXS("../7Li_dp_testL1/Outputs/xs/g1_xs_lDependent.dat")};
    auto* g_merged_1st {MergeGraphs(g_L1_xs_1st, g_sil_xs_1st)};

    Angular::Comparator comp {"Merged g.s", g_merged};
    comp.Add("ADWA", "../7Li_dp/Inputs/gs_ADWA/fort.202");
    comp.Add("ADWA DWUCK", "../FranckCalculations/DWUCK4_personal_calculations/good_calculations/gs/xs_dw4.n");
    comp.Fit();
    comp.Draw("", true);
    Angular::Comparator compSil {"Sil g.s", g_sil_xs};
    compSil.Add("DA1pcorr", "./Inputs/gs_DA1pcorr_Delaroche/21.g0");
    compSil.Add("ADWA", "../7Li_dp/Inputs/gs_ADWA/fort.202");
    compSil.Add("ADWA DWUCK", "../FranckCalculations/DWUCK4_personal_calculations/good_calculations/gs/xs_dw4.n");
    compSil.Fit();
    compSil.Draw("", true);
    Angular::Comparator compL1 {"L1 g.s", g_L1_xs};
    compL1.Add("DA1pcorr", "./Inputs/gs_DA1pcorr_Delaroche/21.g0");
    compL1.Add("ADWA", "../7Li_dp/Inputs/gs_ADWA/fort.202");
    compL1.Add("ADWA DWUCK", "../FranckCalculations/DWUCK4_personal_calculations/good_calculations/gs/xs_dw4.n");
    compL1.Fit();
    compL1.Draw("", true);

    Angular::Comparator comp1 {"Merged 1st Ex", g_merged_1st};
    comp1.Add("DA1pcorr", "./Inputs/g1_DA1pcorr_Delaroche/21.g1");
    comp1.Add("ADWA", "../7Li_dp/Inputs/gs_ADWA/fort.203");
    comp1.Add("ADWA DWUCK", "../FranckCalculations/DWUCK4_personal_calculations/good_calculations/g1/xs_dw4.n");
    comp1.Fit();
    comp1.Draw("", true);
    Angular::Comparator compSil1 {"Sil 1st Ex", g_sil_xs_1st};
    compSil1.Add("DA1pcorr", "./Inputs/g1_DA1pcorr_Delaroche/21.g1");
    compSil1.Add("ADWA", "../7Li_dp/Inputs/gs_ADWA/fort.203");
    compSil1.Add("ADWA DWUCK", "../FranckCalculations/DWUCK4_personal_calculations/good_calculations/g1/xs_dw4.n");
    compSil1.Fit();
    compSil1.Draw("", true);
    Angular::Comparator compL11 {"L1 1st Ex", g_L1_xs_1st};
    compL11.Add("DA1pcorr", "./Inputs/g1_DA1pcorr_Delaroche/21.g1");
    compL11.Add("ADWA", "../7Li_dp/Inputs/gs_ADWA/fort.203");
    compL11.Add("ADWA DWUCK", "../FranckCalculations/DWUCK4_personal_calculations/good_calculations/g1/xs_dw4.n");
    compL11.Fit();
    compL11.Draw("", true);

    // ---------- v0 ... v8 (estados no ligados) ----------
    for(int i = 0; i <= 8; ++i)
        ProcessUnbound("v" + std::to_string(i));
}