#include "ActDataManager.h"
#include "ActInputParser.h"
#include "ActLine.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActTPCData.h"
#include "ActTypes.h"

#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"

#include <iostream>
#include <vector>

typedef ROOT::RVecF Vector;

void Get()
{
    // Read data
    ROOT::EnableImplicitMT();
    std::string beam {"11Li"};
    ActRoot::DataManager datman {"../../configs/data_" + beam + ".conf", ActRoot::ModeType::EFilter};
    std::string moment {"post_preMeasure"};
    if(beam == "7Li")
        datman.SetRuns(66, 69); // 7Li runs after change in L0 trigger in run 67
    else if(beam == "11Li")
    {
        if(moment == "pre")
            datman.SetRuns(19, 65); // 11Li runs pre 7Li
        else if(moment == "post_preMeasure")
            datman.SetRuns(95, 107); // 11Li runs post 7Li
        else if(moment == "post_postMeasure")
            datman.SetRuns(108, 117); // 11Li runs post 7Li
        // else if(moment == "test")
        //     datman.SetRuns(107, 108); // 11Li runs post 7Li
        // datman.SetRuns(19, 65); // 11Li runs pre 7Li
        // datman.SetRuns(95, 122); // 11Li runs post 7Li
        // else if(moment == "all") // in this case do not set any run
    }

    auto chain {datman.GetChain()};
    auto chain2 {datman.GetChain(ActRoot::ModeType::EMerge)};
    auto chain3 {datman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(chain2.get());
    chain->AddFriend(chain3.get());
    chain->SetBranchStatus("fClusters.fVoxels", false);
    ROOT::RDataFrame df {*chain};

    // Read conversion factors
    float padSide {2}; // mm
    ActRoot::InputParser parser {"../../configs/detector.conf"};
    auto merger {parser.GetBlock("Merger")};
    float driftFactor {static_cast<float>(merger->GetDouble("DriftFactor"))};
    std::cout << "-> DriftFactor : " << driftFactor << '\n';

    // Filter by GATCONF and only one BL in cluster vector
    auto def {df.Filter([](ActRoot::ModularData& m) { return m.Get("GATCONF") == 64; }, {"ModularData"})
                  .Filter("fClusters.size() == 1")
                  .Filter("fClusters.fIsBeamLike.front() == true")
                  .Filter(
                      [](const Vector& first, const Vector& second)
                      {
                          bool condA {first.front() <= 0};
                          bool condB {second.front() >= 127};
                          return condA && condB;
                      },
                      {"fClusters.fXRange.first", "fClusters.fXRange.second"})
                  .Define("Line",
                          [&](ActRoot::TPCData& data, ActRoot::MergerData& m)
                          {
                              auto& voxels {data.fClusters.front().GetVoxels()};
                              ActRoot::Line line;
                              line.FitVoxels(voxels, true, true, true);
                              // correction for change in L0 trigger window
                              if(m.fRun < 67)
                              {
                                  auto p = line.GetPoint();
                                  p.SetZ(p.Z() + 62.5 / 4.); // Offset of 5 microsecons == 62.5 tb
                                  line.SetPoint(p);
                              }
                              // conversion to physical units
                              line.Scale(padSide, driftFactor);
                              return line;
                          },
                          {"TPCData", "MergerData"})
                  .Define("AtBegin", [](const ActRoot::Line& l) { return l.MoveToX(0); }, {"Line"})
                  .Define("AtEnd", [](const ActRoot::Line& l) { return l.MoveToX(256); }, {"Line"})};
    auto count {def.Count()};
    ROOT::RDF::Experimental::AddProgressBar(def);
    def.Snapshot("Emittance_Tree", "./Outputs/emittance" + beam + "_" + moment + ".root",
                 {"AtBegin", "AtEnd", "Line", "fRun", "fEntry"});
    std::cout << "Processed events : " << *count << '\n';

    // ROOT::DisableImplicitMT();

    // Get events into file .dat to debug two blobs on 7Li
    // std::ofstream streamer {"./Outputs/events.dat"};
    // def.Foreach(
    //     [&](ActRoot::MergerData& m)
    //     {
    //         m.Stream(streamer);
    //     },
    //     {"MergerData"});
    // streamer.close();

    // Get events with shift in z
    // auto bad = def.Filter(
    // [](const ActRoot::Line& l)
    // {
    //     return l.GetPoint().Z() < 200;
    // },
    // {"Line"});
    //
    // std::ofstream streamer {"./Outputs/eventsVerticalShift.dat"};
    // bad.Foreach(
    //     [&](ActRoot::MergerData& m)
    //     {
    //         m.Stream(streamer);
    //     },
    //     {"MergerData"});
    // streamer.close();
}