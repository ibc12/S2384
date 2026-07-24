#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"

#include <fstream>

void DebugL1EventsFrontTrigger()
{
    ROOT::EnableImplicitMT();

    // Get file from GetL1EventsFrontTrigger output macro
    auto df {ROOT::RDataFrame("Filtered_Tree", "./Outputs/eventsL1_TriggerF0.root")};

    std::cout << "Number of potential L1 events with front trigger: " << df.Count().GetValue() << std::endl;

    // There are still a lot of events that have charge problems, so let's apply the filtering to those events:
    auto defFiltered = df.Filter(
        [&](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
        {
            auto rp {tpc.fRPs.front()};
            auto rp_x {rp.X()};
            auto rp_y {rp.Y()};
            // Run for all clusters
            int counter {};
            for(auto& cluster : tpc.fClusters)
            {
                auto voxels {cluster.GetRefToVoxels()};
                for(auto& v : voxels)
                {
                    if(v.GetPosition().Y() > rp_y - 5 && v.GetPosition().Y() < rp_y + 5 &&
                       v.GetPosition().X() > rp_x - 10 && v.GetPosition().X() < rp_x + 10) // aprox L1 exclusion zone
                        if(v.GetCharge() > 3000.)
                            counter++;
                }
            }
            if(counter > 1)
                return false;
            return true;
        },
        {"MergerData", "TPCData"});

    // Let's see how many events are left after filtering
    auto nEntriesFiltered = static_cast<int>(defFiltered.Count().GetValue());
    std::cout << "Number of entries after filtering charge around RP: " << nEntriesFiltered << std::endl;

    auto defLight =
        defFiltered
            .Define("lightIdx",
                    [](ActRoot::TPCData& tpc)
                    {
                        auto sizeClusters = tpc.fClusters.size();
                        int lightIdx {0};
                        double maxAngle {-1111};
                        ActRoot::Cluster lightCluster {};
                        ROOT::Math::XYZVectorF beamlikeDir {1, 0, 0};
                        for(int i = 0; i < sizeClusters; ++i)
                        {
                            const auto& cluster = tpc.fClusters[i];
                            double angle =
                                std::abs(TMath::ACos(cluster.GetLine().GetDirection().Unit().Dot(beamlikeDir)));
                            if(angle > maxAngle)
                            {
                                maxAngle = angle;
                                lightCluster = cluster;
                                lightIdx = i;
                            }
                        }
                        return lightIdx;
                    },
                    {"TPCData"})
            .Define("thetaLight",
                    [](ActRoot::TPCData& tpc, int lightIdx)
                    {
                        auto lightCluster = tpc.fClusters[lightIdx];
                        auto lightLine = lightCluster.GetLine();
                        double thetaLight =
                            std::abs(TMath::ACos(lightLine.GetDirection().Unit().Dot(ROOT::Math::XYZVectorF(1, 0, 0))));
                        return thetaLight * TMath::RadToDeg();
                    },
                    {"TPCData", "lightIdx"});

    auto hThetaLight =
        defLight.Histo1D({"hThetaLight", "Theta Light;#theta_{light} (deg);Counts", 100, 0, 180}, "thetaLight");

    auto* c = new TCanvas("c", "c", 800, 600);
    hThetaLight->DrawClone();


    std::ofstream outFile("./Outputs/eventsL1_TriggerF0_debug.dat");
    defFiltered.Foreach([&outFile](ActRoot::MergerData& m) { m.Stream(outFile); }, {"MergerData"});
}