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

struct L1Vars
{
    double fTL {-1};
    double fQtotal {-1};
    double fThetaLight {-1};
};

bool IsInsideActar(ROOT::Math::XYZPointF endPoint, ActRoot::TPCParameters& tpcPars)
{
    int nPadOffset {7};
    bool isInX {endPoint.X() > nPadOffset && endPoint.X() < (tpcPars.GetNPADSX() - nPadOffset)};
    bool isInY {endPoint.Y() > nPadOffset && endPoint.Y() < (tpcPars.GetNPADSY() - nPadOffset)};
    bool isInZ {endPoint.Z() > nPadOffset && endPoint.Z() < (tpcPars.GetNPADSZ() - nPadOffset)};
    return (isInX && isInY && isInZ);
}

void checkL1Events()
{
    ActRoot::DataManager dataman {"../configs/data.conf", ActRoot::ModeType::EMerge};
    auto chain {dataman.GetChain()};
    // Add friends if necessary
    auto friend1 {dataman.GetChain(ActRoot::ModeType::EFilter)};
    chain->AddFriend(friend1.get());
    auto friend2 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(friend2.get());

    // ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {*chain};

    auto tpcPars {ActRoot::TPCParameters("Actar")};
    // Count the number of events
    auto dfFilterL1 = df.Filter(
        [&tpcPars](ActRoot::MergerData& d, ActRoot::ModularData& mod, ActRoot::TPCData& tpc)
        {
            if(mod.Get("GATCONF") == 8)
            {
                if(d.fLightIdx >= tpc.fClusters.size())
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
                return IsInsideActar(endPoint, tpcPars);
            }
            else
                return false;
        },
        {"MergerData", "ModularData", "TPCData"});

    // std::ofstream outFile("./Outputs/L1_events.dat");
    // dfFilterL1.Foreach([&](ActRoot::MergerData &m)
    //                    { m.Stream(outFile); }, {"MergerData"});
    // outFile.close();

    auto count = dfFilterL1.Count();
    std::cout << "Number of L1 events: " << count.GetValue() << std::endl;

    // Define new variables
    auto dfAll = dfFilterL1.Define("L1Vars",
                                   [](ActRoot::MergerData& d, ActRoot::TPCData& tpc)
                                   {
                                       auto& light = tpc.fClusters[d.fLightIdx];
                                       auto& voxels = light.GetRefToVoxels();

                                       auto begin = voxels.front().GetPosition();
                                       auto end = voxels.back().GetPosition();

                                       auto distance {(begin - end).R()};
                                       auto theta3Lab {d.fThetaLight};
                                       double qTotal {};
                                       for(int i = 0; i < voxels.size(); i++)
                                       {
                                           qTotal += voxels[i].GetCharge();
                                       }
                                       L1Vars vars {.fTL = distance, .fQtotal = qTotal, .fThetaLight = theta3Lab};

                                       return vars;
                                   },
                                   {"MergerData", "TPCData"});

    // If cuts are present, apply them
    ActRoot::CutsManager<std::string> cuts;
    // Gas PID
    cuts.ReadCut("L1_p", "./Outputs/p_events_L1.root");

    // Apply PID and save in file
    auto gated {dfAll.Filter(
    [&](const L1Vars& vars) {
        if (cuts.GetCut("L1_p"))
            return cuts.IsInside("L1_p", vars.fTL, vars.fQtotal);
        else
            return false;
    },
    {"L1Vars"})};
    std::ofstream outFile("./Outputs/L1_events_run_22.dat");
    dfAll.Foreach([&](ActRoot::MergerData &m)
                       { m.Stream(outFile); }, {"MergerData"});
    outFile.close();

    auto hitsLengthCharge {dfAll.Define("x", "L1Vars.fTL")
                               .Define("y", "L1Vars.fQtotal")
                               .Histo2D({"hitsLengthCharge", "Length vs Charge;Track Length [mm];Total Charge [u.a]",
                                         120, 0, 100, 1000, 0, 300000},
                                        "x", "y")};
    auto histLengthAngle {dfAll.Define("x", "L1Vars.fTL")
                              .Define("y", "L1Vars.fThetaLight")
                              .Histo2D({"hitsLengthAngle", "Length vs Angle;Track Length [mm];Theta Light [deg]", 120,
                                        0, 100, 100, 0, 180},
                                       "x", "y")};

    // Create a canvas to visualize the results
    TCanvas* c1 = new TCanvas("c1", "L1 Events Count", 800, 600);
    c1->Divide(2, 1);
    c1->cd(1);
    hitsLengthCharge->DrawClone("COLZ");
    c1->cd(2);
    histLengthAngle->DrawClone("COLZ");
}