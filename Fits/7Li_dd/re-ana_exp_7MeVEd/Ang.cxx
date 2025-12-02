#include "ActKinematics.h"
#include "ActMergerData.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TColor.h"
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

#include <string>
#include <vector>

#include "../../Histos.h"

// void Angular::Comparator::Flip(const std::string& name)
// {
//     auto graphs = this->GetTheoGraphs();
//     // verify graph exists
//     if(graphs.find(name) == graphs.end())
//         throw std::runtime_error("Comparator::Flip - Graph '" + name + "' not found");
//
//     // Get original graph and clone to avoid double-delete
//     auto gOrig = graphs.at(name);
//     auto gFlip = (TGraphErrors*) gOrig->Clone((name + "_flipped").c_str());
//
//     // Aply transformation θ' = 180 - θ
//     for(int i = 0; i < gFlip->GetN(); ++i)
//     {
//         double x, y;
//         gFlip->GetPoint(i, x, y);
//         gFlip->SetPoint(i, 180.0 - x, y);
//     }
//     // Sort from smallest to largest to be able to plot and replace it
//     gFlip->Sort();
//     this->Replace(name, gFlip);
// }

void Ang(bool isLab = false)
{

    if(isLab)
        Angular::ToggleIsLab();

    // Experimental data
    TGraphErrors* gExp_7MeV_Ed {new TGraphErrors("./Inputs/14-7MeVEd_paper_DA1p.dat", "%lg %lg")};

    // Plot
    Angular::Comparator comp {"g.s", gExp_7MeV_Ed};
    // comp.Add("DA1p", "./Inputs/gsDA1p_14-7Ed/fort.201");
    comp.Add("OMP paper", "../Inputs/gsDA1p_corr/fort.201");
    //auto graphs = comp.GetTheoGraphs();
    // auto gOrig = graphs.at("OMP paper");
    //
    //// Clonar
    // auto gFlip = (TGraphErrors*)gOrig->Clone("paper_flipped");
    //
    //// Flip angular
    // for(int i = 0; i < gFlip->GetN(); ++i)
    //{
    //     double x, y;
    //     gFlip->GetPoint(i, x, y);
    //     gFlip->SetPoint(i, 180.0 - x, y);
    // }
    //
    //// Ordenar para que quede 0 → 180 creciente
    // gFlip->Sort();
    //
    //// Reemplazar
    // comp.Replace("OMP paper", gFlip);
    comp.Fit();
    comp.Draw("", true);
}
