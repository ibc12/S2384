#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActParticle.h"
#include "ActSRIM.h"
#include "ActSilData.h"
#include "ActSilMatrix.h"
#include "ActSilSpecs.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TH2.h"
#include "TH2D.h"
#include "TString.h"

#include <fstream>
#include <iostream>
#include <map>
#include <string>

#include "./../../PrettyStyle.C"

void GetDiffusionParameters()
{
    // Get diffusion parameters from simulation
    // PrettyStyle(false);
    std::string dataconf {"../../configs/data.conf"};

    // Read data
    ActRoot::DataManager dataman {dataconf, ActRoot::ModeType::EMerge};
    dataman.SetRuns(64, 67);
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    auto chain3 {dataman.GetChain(ActRoot::ModeType::EFilter)};
    chain->AddFriend(chain2.get());
    chain->AddFriend(chain3.get());

    // RDataFrame
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {*chain};

    // Get the L1 events
    auto lambdaIsL1 {[](ActRoot::MergerData& mer, ActRoot::ModularData& mod)
                     { return (mod.Get("GATCONF") == 8) && (mer.fLightIdx != -1); }};
    auto dfL1 = df.Filter(lambdaIsL1, {"MergerData", "ModularData"});

    // Go for events with z almost continuous
    // Filter low charge deposition (high energy light particles, RANSAC needed)
    auto dfFilterL1 = dfL1.Filter(
        [](ActRoot::MergerData& mer, ActRoot::TPCData& tpc)
        {
            auto thetaLight = mer.fThetaLight;
            auto phiLight = std::abs(mer.fPhiLight); // aboid -90 deg ambiguity
            if (phiLight < 100.0 && phiLight > 80.0 && mer.fLightIdx != -1)
                return true;
            else if (phiLight > 260.0 && phiLight < 280.0 && mer.fLightIdx != -1)
                return true;
            else
                return false;
        },
        {"MergerData", "TPCData"}).Filter(
        [](ActRoot::MergerData& mer)
        {
            if (mer.fLight.fQave > 300.0)
                return true;
            else
                return false;
        },
        {"MergerData"});

    // Look at the events to ensure they are good candidates
    std::ofstream outFile {"./Outputs/L1events_run_64_67.dat"};
    dfFilterL1.Foreach([&outFile](const ActRoot::MergerData& mer) { mer.Stream(outFile); }, {"MergerData"});
    outFile.close();
}