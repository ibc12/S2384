#include "ActCutsManager.h"
#include "ActMergerData.h"

#include <ROOT/RDataFrame.hxx>

#include <fstream>

void GetEventsInsideCut()
{
    ROOT::DisableImplicitMT();
    ROOT::RDataFrame df {"Final_Tree", "../PostAnalysis/Outputs/tree_ex_F_11Li_d_p_filtered.root"};

    // Get cut in kinematics
    ActRoot::CutsManager<std::string> cuts;
    // kienmatic cut
    cuts.ReadCut("events", TString::Format("../PostAnalysis/Cuts/pid_p_l1_11Li.root").Data());

    // Filter dataframe with cut
    auto gated {df.Filter([&cuts](ActRoot::MergerData& mer, const double EVertex)
                          { return cuts.IsInside("events", mer.fLight.fRawTL, mer.fLight.fQtotal); }, {"MergerData", "EVertex"})};

    // Save the filtered dataframe events
    std::ofstream out("./Outputs/p_events_L1_12Li.dat");
    gated.Foreach([&](ActRoot::MergerData& m) { m.Stream(out); }, {"MergerData"});
    out.close();
}