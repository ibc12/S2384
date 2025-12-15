#include "ActCutsManager.h"
#include "ActMergerData.h"

#include <ROOT/RDataFrame.hxx>

#include <fstream>

void GetEventsInsideCut()
{
    ROOT::DisableImplicitMT();
    ROOT::RDataFrame df {"Final_Tree", "../PostAnalysis/Outputs/tree_ex_11Li_d_p_filtered.root"};

    // Get cut in kinematics
    ActRoot::CutsManager<std::string> cuts;
    // kienmatic cut
    cuts.ReadCut("events", TString::Format("../PostAnalysis/eventsWithStructure_12Li.root").Data());

    // Filter dataframe with cut
    auto gated {df.Filter([&cuts](ActRoot::MergerData& mer, const double EVertex)
                          { return cuts.IsInside("events", mer.fThetaLight, EVertex); }, {"MergerData", "EVertex"})};

    // Save the filtered dataframe events
    std::ofstream out("./Outputs/eventsWithStructureInKinematics_12Li.dat");
    gated.Foreach([&](ActRoot::MergerData& m) { m.Stream(out); }, {"MergerData"});
    out.close();
}