#include "ActMergerData.h"

#include "ROOT/RDataFrame.hxx"


void CheckEventsBadFiltering()
{
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"PreProcessed_Tree", "../PostAnalysis/Outputs/tree_preprocess_11Li.root"};

    auto dfFilter = df.Filter([](ActRoot::MergerData& m) { return m.fEntry == 58803 && m.fRun == 21; }, {"MergerData"});

    dfFilter.Foreach(
        [](ActRoot::MergerData& m)
        {
            auto rp = m.fRP.X();
            auto theta = m.fThetaHeavy;
            auto TL = m.fHeavy.fTL;

            std::cout << "rp at: " << rp << ", theta: " << theta << ", TL: " << TL << std::endl;
            auto finalPosition = rp + TL * TMath::Cos(theta * TMath::DegToRad());
            std::cout << "Final position: " << finalPosition << std::endl;
        },
        {"MergerData"});
}