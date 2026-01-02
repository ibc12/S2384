#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActSilData.h"
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

#include "../../PostAnalysis/HistConfig.h"

void FrontSilCalibrationCheck()
{
    // Let's check the channel hitted for the front silicons. Geater part of hits seem to have very big energy. Indeed,
    // the energy is greater that the punchthough for what I think is Z=1 particles. Let's start by the t selection on
    // 11Li beam

    //  Read data
    auto filename {TString::Format("./Inputs/tree_ex_11Li_d_t_filtered.root")};
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"Final_Tree", filename};

    auto def = df.Define("SilE", [&](const ActRoot::SilData& s)
                     { return s.fSiE.at("f0").front(); }, {"SilData"});

    // Plot all silicons together (similar gain)
    auto hSilE = def.Histo1D({"SilE", "SilE;SilE [channel]", 17000, 0, 17000}, "SilE");

    auto* c {new TCanvas("c", "c")};
    c->cd();
    hSilE->DrawClone();
}