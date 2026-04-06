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

#include <cstdint>
#include <fstream>
#include <map>
#include <mutex>
#include <string>
#include <unordered_set>
#include <vector>

void checkPIDdifferentRuns()
{
    // Thomas told me that maybe there was a gain change in the middle of the experiment, so let's plot the PID in
    // intervals of 5 runs to see if there is any change


    // 1. Get the data
    std::string beam {"11Li"};
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

    // RDataFrame
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame dforigin {*chain};

    int interval = 5;
    int run_min = 15;
    int run_max = 125;

    // 2. Create the histograms and fill them
    std::map<int, ROOT::TThreadedObject<TH2D>> hsPID;
    for(int run = run_min; run <= run_max; run += interval)
    {
        int run_5 = run / interval;
        hsPID.emplace(run_5, TH2D(TString::Format("hPID_run%d", run_5),
                                  TString::Format("PID for run %d to %d;E_{Sil} [MeV];#Delta E_{gas} [arb. units]", run, run + interval - 1),
                                  450, 0, 70, 600, 0, 3000));
    }

    auto lambdaPID = [](ActRoot::MergerData& m)
    {
        if(m.fLight.GetNLayers() == 1)
            return m.fLight.GetLayer(0) == "f0";
        else
            return false;
    };
    auto dfPID = dforigin.Filter(lambdaPID, {"MergerData"});
    dfPID.Foreach([&](ActRoot::MergerData& m)
                  { hsPID[m.fRun / interval].Get()->Fill(m.fLight.fEs[0], m.fLight.fQave); }, {"MergerData"});

    // 3. Plot them
    auto* cPID = new TCanvas("cPID", "PID for different runs", 1200, 800);
    cPID->Divide(5, 5);
    int pad = 1;
    for(int run = run_min; run <= run_max; run += interval)
    {
        cPID->cd(pad);

        auto h = hsPID[run / interval].Merge();
        if(h)
            h->DrawClone("colz");

        pad++;
    }
}