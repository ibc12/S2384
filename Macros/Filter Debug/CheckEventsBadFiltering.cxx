#include "ActDataManager.h"
#include "ActMergerData.h"

#include "ROOT/RDataFrame.hxx"


void CheckEventsBadFiltering()
{
    ROOT::EnableImplicitMT();

    ActRoot::DataManager dataManager {};
    dataManager.ReadDataFile("../../configs/data.conf");

    // For L1 enough with 4 runs (64, 67). For lat sils put at least 10
    dataManager.SetRuns(20, 30);

    // Get df for the runs
    auto chain {dataManager.GetChain(ActRoot::ModeType::EReadSilMod)};
    auto chain1 {dataManager.GetChain(ActRoot::ModeType::EMerge)};
    chain->AddFriend(chain1.get());
    auto dfAll {ROOT::RDataFrame(*chain)};

    dfAll.Foreach(
        [](ActRoot::MergerData& m)
        {
            if((m.fRun == 30 && m.fEntry == 123312) || (m.fRun == 23 && m.fEntry == 28439))
            {
                auto rp = m.fRP.X();
                auto theta = m.fThetaHeavy;
                auto TL = m.fHeavy.fTL;

                auto finalPosition = rp + TL * TMath::Cos(theta * TMath::DegToRad());
                // if(finalPosition < 240)
                std::cout << "All:   Run: " << m.fRun << " Event: " << m.fEntry << " x_end: " << finalPosition
                          << " rp.x: " << m.fRP.X() << " TL: " << TL << " thetaMerger: " << m.fThetaHeavy
                          << " thetaLight: " << m.fThetaLight << std::endl;
            }
        },
        {"MergerData"});

    ROOT::RDataFrame df {"PreProcessed_Tree", "../../PostAnalysis/Outputs/tree_preprocess_11Li.root"};

    // auto dfFilter = df.Filter([](ActRoot::MergerData& m) { return m.fEntry == 58803 && m.fRun == 21; },
    // {"MergerData"});

    df.Foreach(
        [](ActRoot::MergerData& m)
        {
            if((m.fRun == 30 && m.fEntry == 123312) || (m.fRun == 23 && m.fEntry == 28439))
            {
                auto rp = m.fRP.X();
                auto theta = m.fThetaHeavy;
                auto TL = m.fHeavy.fTL;

                auto finalPosition = rp + TL * TMath::Cos(theta * TMath::DegToRad());

                std::cout << "Run: " << m.fRun << " Event: " << m.fEntry << " x_end: " << finalPosition
                          << " rp.x: " << m.fRP.X() << " TL: " << TL << " thetaMerger: " << m.fThetaHeavy
                          << " thetaLight: " << m.fThetaLight << std::endl;
            }
        },
        {"MergerData"});
}