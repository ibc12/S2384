#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActSilData.h"
#include "ActTPCData.h"
#include "ActTPCParameters.h"
#include "ActCutsManager.h"

#include <ROOT/RDataFrame.hxx>

#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"

#include <fstream>
#include <iostream>


bool IsInsideActar(ROOT::Math::XYZPointF endPoint, ActRoot::TPCParameters& tpcPars)
{
    int nPadOffset {2};
    bool isInX {endPoint.X() > nPadOffset && endPoint.X() < (tpcPars.GetNPADSX() - nPadOffset)};
    bool isInY {endPoint.Y() > nPadOffset && endPoint.Y() < (tpcPars.GetNPADSY() - nPadOffset)};
    bool isInZ {endPoint.Z() > nPadOffset && endPoint.Z() < (tpcPars.GetNPADSZ() - nPadOffset)};
    return (isInX && isInY && isInZ);
}

bool IsAtEndZ(ROOT::Math::XYZPointF endPoint, ActRoot::TPCParameters& tpcPars)
{
    return (endPoint.Z() >= (tpcPars.GetNPADSZ() - 10));
}

void checkIfBrokenTracksOnSilEvents()
{
    ActRoot::DataManager dataman {"../configs/data.conf", ActRoot::ModeType::EMerge};
    auto chain {dataman.GetChain()};
    // Add friends if necessary
    auto friend1 {dataman.GetChain(ActRoot::ModeType::EFilter)};
    chain->AddFriend(friend1.get());
    auto friend2 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(friend2.get());

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {*chain};

    df.Describe().Print();

    auto tpcPars {ActRoot::TPCParameters("Actar")};
    // Count the number of events
    auto dfFilter = df.Filter(
        [&tpcPars](ActRoot::MergerData& d, ActRoot::ModularData& mod, ActRoot::TPCData& tpc)
        {
            if(mod.Get("GATCONF") == 1)
            {
                if(d.fLightIdx == -1)
                    return false;

                auto& light = tpc.fClusters[d.fLightIdx];
                auto direction = light.GetLine().GetDirection();
                light.SortAlongDir(direction);
                auto& voxels = light.GetRefToVoxels();
                if(voxels.size() < 2)
                    return false;

                auto begin = voxels.front().GetPosition();
                auto end = voxels.back().GetPosition();

                if(direction.R() == 0)
                    return false;

                auto trackLength = (begin - end).R();
                auto endPoint = d.fRP + trackLength * direction.Unit();
                return IsAtEndZ(endPoint, tpcPars);
            }
            else
                return false;
        },
        {"MergerData", "ModularData", "TPCData"});

    auto dfSave = df.Filter(
        [](ActRoot::MergerData& d, ActRoot::ModularData& mod)
        {
            if(mod.Get("GATCONF") == 1 && d.fPhiLight == -1)
                return true;
            else
                return false;
        },
        {"MergerData", "ModularData"});
    // std::ofstream outFile("./Outputs/Phi_Minus1_LeftSil_events.dat");
    // dfSave.Foreach([&](ActRoot::MergerData &m)
    //                    { m.Stream(outFile); }, {"MergerData"});
    // outFile.close();

    auto count = dfFilter.Count();
    std::cout << "Number of sil left events not reaching out: " << count.GetValue() << std::endl;

    
    auto histPhi {df.Histo1D({"histPhi", "Phi;#phi [rad]", 100, -TMath::Pi(), TMath::Pi()}, "fPhiLight")};

    // Create a canvas to visualize the results
    TCanvas* c1 = new TCanvas("c1", "L1 Events Count", 800, 600);
    c1->Divide(2, 1);
    c1->cd(1);
    histPhi->DrawClone("colz");
}