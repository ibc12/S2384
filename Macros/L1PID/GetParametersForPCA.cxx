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

void GetParametersForPCA()
{
    PrettyStyle(true, true);
    // Get all data for 7Li (there is no triton there, so easier to see deuterium)
    std::string beam {"11Li"};
    std::string target {"d"};
    std::string light {"d"};

    std::string dataconf {};
    if(beam == "11Li")
        dataconf = "../../configs/data_11Li.conf";
    else if(beam == "7Li")
        dataconf = "../../configs/data_7Li.conf";
    else
        throw std::runtime_error("Beam cannot differ from 11Li or 7Li");

    // Read data
    ActRoot::DataManager dataman {dataconf, ActRoot::ModeType::EMerge};
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(chain2.get());

    // RDataFrame
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame dforigin {*chain};

    // Filter GATCONF == L1
    auto df =
        dforigin.Filter([](ActRoot::ModularData& mod, ActRoot::MergerData& mer)
                        { return mod.Get("GATCONF") == 8 && (mer.fLightIdx != -1); }, {"ModularData", "MergerData"});

    df.Describe().Print();
    // Snapshot filtered dataframe Qave ThetaLight Qtot TrackLength
    df.Snapshot("L1PCA", "./Outputs/L1PCA_11Li.root",
                {"fQave", "fLight.fQtotal", "fThetaLight", "fLight.fTL", "fPhiLight"});
}