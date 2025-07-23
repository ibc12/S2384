#include "ActDataManager.h"
#include "ActTPCData.h"

#include <ROOT/RDataFrame.hxx>

#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"

void checkBeamEmitance()
{
    ActRoot::DataManager dataman {"../configs/data.conf", ActRoot::ModeType::EFilter};
    auto chain {dataman.GetChain()};

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {*chain};
    df.Describe().Print();

    // Get data from filtered data
    auto file {new TFile("../RootFiles/Data/Filter_Run_0004.root")};
    auto df {new ROOT::RDataFrame("FilterData", file)};

    // Filter only one cluster events
    auto dfFilter = df.Filter([](ActRoot::TPCData& d) { return d.fClusters.size() == 1; }, {"TPCData"});

    // Define Y and Z variables in window and end of chamber
    dfFilter = dfFilter.Define("yWindow", [](ActRoot::TPCData& d) { return d.fClusters[0].GetLine().MoveToX(0).Y(); })
                   .Define("zWindow", [](ActRoot::TPCData& d) { return d.fClusters[0].GetLine().MoveToX(0).Z(); })
                   .Define("yEnd", [](ActRoot::TPCData& d) { return d.fClusters[0].GetLine().MoveToX(128).Y(); })
                   .Define("zEnd", [](ActRoot::TPCData& d) { return d.fClusters[0].GetLine().MoveToX(128).Z(); })
                   .Define("yAngle",
                           [](ActRoot::TPCData& d)
                           {
                               auto directionModule {d.fClusters[0].GetLine().GetDirection().R()};
                               return TMath::ACos(d.fClusters[0].GetLine().GetDirection().Y() / directionModule) *
                                      TMath::RadToDeg();
                           })
                   .Define("zAngle",
                           [](ActRoot::TPCData& d)
                           {
                               auto directionModule {d.fClusters[0].GetLine().GetDirection().R()};
                               return TMath::ACos(d.fClusters[0].GetLine().GetDirection().Z() / directionModule) *
                                      TMath::RadToDeg();
                           });

    // Create the histogram
    auto hWindowPosition = dfFilter.Histo2D(
        {"hWindowPosition", "Window Position;Y [cm];Z [cm]", 100, 0, 128, 100, 0, 128}, "yWindow", "zWindow");
    auto hEndPosition =
        dfFilter.Histo2D({"hEndPosition", "End Position;Y [cm];Z [cm]", 100, 0, 128, 100, 0, 128}, "yEnd", "zEnd");
    auto hAngleY = dfFilter.Histo1D({"hAngleY", "Y Angle;Angle [deg];Counts", 100, -10, 10}, "yAngle");
    auto hAngleZ = dfFilter.Histo1D({"hAngleZ", "Z Angle;Angle [deg];Counts", 100, -10, 10}, "zAngle");

    // Create canvas and draw histograms
    TCanvas* c1 = new TCanvas("c1", "Beam Emitance", 1200, 800);
    c1->Divide(2, 2);
    c1->cd(1);
    hWindowPosition->Draw("COLZ");
    c1->cd(2);
    hEndPosition->Draw("COLZ");
    c1->cd(3);
    hAngleY->Draw();
    c1->cd(4);
    hAngleZ->Draw();
}