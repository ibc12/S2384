#include "ActDataManager.h"
#include "ActSilData.h"
#include "ActTPCData.h"
#include "ActModularData.h"
#include "ActMergerData.h"
#include "ActTPCParameters.h"

#include <ROOT/RDataFrame.hxx>

#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"

bool IsInsideActar(ROOT::Math::XYZPointF endPoint, ActRoot::TPCParameters &tpcPars)
{
    bool isInX{endPoint.X() > 0 && endPoint.X() < tpcPars.GetNPADSX()};
    bool isInY{endPoint.Y() > 0 && endPoint.Y() < tpcPars.GetNPADSY()};
    bool isInZ{endPoint.Z() > 0 && endPoint.Z() < tpcPars.GetNPADSZ()};
    return (isInX && isInY && isInZ);
}

void checkL1Events()
{
    ActRoot::DataManager dataman{"../configs/data.conf", ActRoot::ModeType::EMerge};
    auto chain{dataman.GetChain()};
    // Add friends if necessary
    auto friend1{dataman.GetChain(ActRoot::ModeType::EFilter)};
    chain->AddFriend(friend1.get());

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df{*chain};

    // Filter for L1 events
    // auto gated{df.Filter([](ActRoot::ModularData &d)
    //                      { return d.Get("GATCONF") == 8; }, {"ModularData"})};

    auto tpcPars{ActRoot::TPCParameters("Actar")};
    // Count the number of events
    auto dfFilterL1 = df.Filter(
        [&tpcPars](ActRoot::MergerData &d, ActRoot::TPCData &tpc)
        {
            if (d.fLightIdx >= tpc.fClusters.size())
                return false;

            auto &light = tpc.fClusters[d.fLightIdx];
            auto &voxels = light.GetRefToVoxels();
            if (voxels.size() < 2)
                return false;

            auto begin = voxels.front().GetPosition();
            auto end = voxels.back().GetPosition();

            auto direction = light.GetLine().GetDirection();
            if (direction.R() == 0)
                return false;

            auto trackLength = (begin - end).R();
            auto endPoint = d.fRP + trackLength * direction.Unit();
            return IsInsideActar(endPoint, tpcPars);
        },
        {"MergerData", "TPCData"});

    

    auto count = dfFilterL1.Count();
    std::cout << "Number of L1 events: " << count.GetValue() << std::endl;

    // Create a canvas to visualize the results
    TCanvas *c1 = new TCanvas("c1", "L1 Events Count", 800, 600);
    c1->cd();
}