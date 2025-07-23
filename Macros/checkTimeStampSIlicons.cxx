
#include "TTree.h"
#include "ActSilData.h"
#include <ROOT/RDataFrame.hxx>
#include "TCanvas.h"
#include "ActModularData.h"

#include "TFile.h"

#include <iostream>
#include <vector>
#include <algorithm>

void checkTimeStampSIlicons()
{
    // Get data from file in TTree
    auto fileSilicon{new TFile("../RootFiles/Data/Data_Run_0035.root")};
    auto fileRaw{new TFile("../RootFiles/Raw/Tree_Run_0035_Merged.root")};

    auto treeSilicon{fileSilicon->Get<TTree>("VXITree")};
    auto treeRaw{fileRaw->Get<TTree>("ACTAR_TTree")};

    treeSilicon->AddFriend("ACTAR_TTree", fileRaw);

    ROOT::RDataFrame dfSilicon(*treeSilicon);

    auto df = dfSilicon.Define("time", [&](ActRoot::ModularData &md)
    {
        auto timeh_up {md.Get("CTR_TIMEH_UP")};
        auto timeh {md.Get("CTR_TIMEH")};
        auto timeml_up {md.Get("CTR_TIMEML_UP")};
        auto timeml {md.Get("CTR_TIMEML")};

        return timeml + timeml_up * 2e16 + timeh * 2e32 + timeh_up * 2e48;
    }, {"ModularData"});

    auto minTime {*df.Min("time")};

    df = df.Define("timeCorrected", [&](double &time)
    {
        // Correct the time
        return (time - minTime) ;
    }, {"time"});

    df = df.Define("Energy", [&](ActRoot::SilData &d)
                               {
                                std::string layer {"l0"};
        if (d.fSiN.count(layer))
        {
            auto& ns {d.fSiN[layer]};
            auto it{std::find(ns.begin(), ns.end(), 0)};
            if (it != ns.end())
            {
                auto index = std::distance(ns.begin(), it);
                return d.fSiE[layer][index];
            }}

            return -1.f; }, {"SilData"});

    // Create the th2 model
    const ROOT::RDF::TH2DModel ESil_TimeStamp {"hl0_0", "l0_0;Nevent; E_{Sil} [channel]", 600, 0, 600e7, 1000, 0,
                                5000};
    auto histo {df.Histo2D(ESil_TimeStamp, "timeCorrected", "Energy")};
    auto hTime {df.Histo1D("timeCorrected")};
    auto hml {df.Define("ml", "fLeaves[\"CTR_TIMEML\"]").Histo1D("ml")};

    auto c {new TCanvas("c", "Silicon Energy vs time")};
    histo->DrawClone("colz");

    auto cTime {new TCanvas("cTime", "Silicon Time")};
    cTime->DivideSquare(4);
    cTime->cd(1);
    hTime->DrawClone("colz");
    cTime->cd(2);
    hml->DrawClone();
}