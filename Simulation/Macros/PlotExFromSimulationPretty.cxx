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

#include "../../PrettyStyle.C"


void PlotExFromSimulationPretty()
{
    PrettyStyle(false, false);

    std::string fileName {"../Outputs/7Li/Decay/2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_decay_democratic4Body.root"};

    ROOT::EnableImplicitMT();

    ROOT::RDataFrame df("SimulationTTree", fileName.c_str());

    // Plot branch Eex
    auto hEx {df.Histo1D({"hEx", "Ex;E_{x} [MeV];Counts", 300, 0, 12}, "Eex")};

    TCanvas* c1 {new TCanvas("c1", "c1", 1600, 1600)};

    // Adjust canvas margins
    c1->SetLeftMargin(0.16);
    c1->SetRightMargin(0.05);
    c1->SetBottomMargin(0.12);
    c1->SetTopMargin(0.08);

    // Move the Y-axis title further away from the tick labels
    hEx->GetYaxis()->SetTitleOffset(1.6);

    hEx->DrawClone();

    // Force ROOT to recalculate the layout before saving
    c1->Modified();
    c1->Update();

    c1->SaveAs("../../../Charlas/Ganil Colloque 2026/Figures presentation/7Li_dp_AlfaT_PhaseSpace.png", "PNG");
}