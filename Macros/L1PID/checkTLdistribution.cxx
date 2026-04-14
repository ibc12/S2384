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

#include "../../PrettyStyle.C"


void checkTLdistribution()
{
    // Get df to filter only deuteriums
    // Get all data for 7Li (there is no triton there, so easier to see deuterium)
    std::string beam {"11Li"};
    std::string target {"d"};
    std::string light {"p"};

    auto inFile {TString::Format("../../PostAnalysis/Outputs/tree_preprocess_F_%s.root", beam.c_str())};

    // RDataFrame
    // ROOT::EnableImplicitMT();
    ROOT::RDataFrame dforigin {"PreProcessed_Tree", inFile};

    // Filter GATCONF == L1
    auto def =
        dforigin.Filter([](ActRoot::ModularData& mod, ActRoot::MergerData& mer)
                        { return mod.Get("GATCONF") == 8 && (mer.fLightIdx != -1); }, {"ModularData", "MergerData"});

    ActRoot::CutsManager<std::string> cuts;
    // Save some events
    // Get only deuteriums, and then gate them in theta for only get the elastic, and also in phi
    cuts.ReadCut("deuterium", "./Cuts/d_TLvsQ_11Li.root");
    auto df = def.Filter(
        [&](ActRoot::MergerData& m)
        {
            if(cuts.IsInside("deuterium", m.fLight.fRawTL, m.fLight.fQtotal))
            {
                return true;
            }
            return false;
        },
        {"MergerData"});
    auto dfPhiPositive = df.Filter([](ActRoot::MergerData& m)
                          { return (m.fPhiLight > 80 && m.fPhiLight < 100) && (m.fThetaLight > 85 && m.fThetaLight < 90); },
                          {"MergerData"});
    auto dfPhiNegative = df.Filter([](ActRoot::MergerData& m)
                          { return (m.fPhiLight > -100 && m.fPhiLight < -80) && (m.fThetaLight > 85 && m.fThetaLight < 90); },
                          {"MergerData"}); 

    auto h_TL_PhiPositiveElastic = dfPhiPositive.Histo1D({"h_TL_PhiPositiveElastic", "TL distribution for deuterium elastic events with phi near 90 positive;Raw TL [a.u.];Counts", 240, 0, 120},
                                     "fLight.fRawTL");
    auto h_TL_PhiNegativeElastic = dfPhiNegative.Histo1D({"h_TL_PhiNegativeElastic", "TL distribution for deuterium elastic events with phi near 90 negative;Raw TL [a.u.];Counts", 240, 0, 120},
                                     "fLight.fRawTL");

    auto* c = new TCanvas("TL distribution for deuterium elastic events", "TL distribution for deuterium elastic events", 1200, 600);
    c->Divide(2, 1);
    c->cd(1);
    h_TL_PhiPositiveElastic->DrawClone();
    c->cd(2);
    h_TL_PhiNegativeElastic->DrawClone();

        
}