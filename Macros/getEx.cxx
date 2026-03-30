#include "ActCluster.h"
#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActLine.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActSRIM.h"
#include "ActTPCData.h"
#include "ActTPCParameters.h"
#include "ActVoxel.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"
#include <random>

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMath.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TSpline.h>
#include <TStyle.h>

#include <Math/Point3D.h>
#include <Math/Vector3D.h>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

#include "../PrettyStyle.C"
#include "../PostAnalysis/HistConfig.h"

void getEx()
{
    std::string beam {"7Li"};
    std::string target {"d"};
    std::string light {"d"};
    // Get two .root files, filter to get silicon events
    auto infile {TString::Format("../PostAnalysis/tree_ex_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str())};
    auto infileFiltered {TString::Format("../PostAnalysis/tree_ex_%s_%s_%s_filtered.root", beam.c_str(), target.c_str(), light.c_str())};

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"Final_Tree", infile.Data()};
    ROOT::RDataFrame dfFiltered {"Final_Tree", infileFiltered.Data()};

    auto dfSil = df.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == true; }, {"MergerData"});
    auto dfFilteredSil = dfFiltered.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == true; }, {"MergerData"});


    auto hExSil = dfSil.Histo1D(HistConfig::ExZoom, "Ex");
    hExSil->SetTitle("Ex with silicons");
    auto hExFilteredSil = dfFilteredSil.Histo1D(HistConfig::ExZoom, "Ex");
    hExFilteredSil->SetTitle("Ex with silicons filtered");

    auto* c = new TCanvas("cEx", "cEx", 800, 600);
    c->Divide(2, 1);
    c->cd(1);
    hExSil->DrawClone();
    c->cd(2);
    hExFilteredSil->DrawClone();

}