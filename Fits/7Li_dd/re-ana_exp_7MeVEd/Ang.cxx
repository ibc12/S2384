#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TROOT.h"
#include "TString.h"

#include "AngComparator.h"
#include "AngDifferentialXS.h"
#include "AngFitter.h"
#include "AngGlobals.h"
#include "AngIntervals.h"
#include "FitInterface.h"
#include "Interpolators.h"
#include "PhysExperiment.h"

#include "ActMergerData.h"
#include "ActKinematics.h"

#include <string>
#include <vector>
#include "TStyle.h"
#include "TColor.h"

#include "../../Histos.h"

void Ang(bool isLab = false)
{

    if(isLab)
        Angular::ToggleIsLab();

    // Experimental data
    TGraphErrors* gExp_7MeV_Ed {new TGraphErrors("./Inputs/7MeVEd.dat", "%lg %lg")};

    // Plot
    Angular::Comparator comp {"g.s", gExp_7MeV_Ed};
    comp.Add("DA1p", "./Inputs/da1p/fort.201");
    comp.Fit();
    comp.Draw("", true);
}
