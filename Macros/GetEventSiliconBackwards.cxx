#include "ActMergerData.h"
#include "ActModularData.h"

#include "ROOT/RDataFrame.hxx"

#include <iostream>
#include <string>
#include <vector>

void GetEventSiliconBackwards()
{
    ROOT::RDataFrame df {"Final_Tree", "../PostAnalysis/Outputs/tree_ex_F_11Li_d_p_filtered.root"};

    auto dfFilter = df.Filter(
        [](ActRoot::MergerData& m)
        {
            return m.fLight.IsFilled() && (m.fLight.fLayers.front() == "l0" || m.fLight.fLayers.front() == "r0") &&
                   m.fThetaLight > 100;
        },
        {"MergerData"});

    dfFilter.Foreach([](ActRoot::MergerData& m)
                     { std::cout << "Event: " << m.fEntry << ", Run: " << m.fRun << std::endl; }, {"MergerData"});
}