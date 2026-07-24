#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"

#include <fstream>

void GetL1EventsFrontTrigger()
{
    ROOT::EnableImplicitMT();

    // ROOT::RDataFrame df {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_7Li_d_d_filtered.root"};

    // // Get the runs from the data manager reading the config file
    ActRoot::DataManager dataManager {};
    dataManager.ReadDataFile("../../configs/data.conf");

    // Get L1ValidationZone parameter
    ActRoot::InputParser parser {};
    parser.ReadFile("../../configs/detector.conf");
    auto mergerBlock = parser.GetBlock("Merger");
    auto validationZoneL1 = mergerBlock->GetDouble("L1ExclusionZone");
    auto driftFactor = mergerBlock->GetDouble("DriftFactor");

    // Per-run zDrift cuts, same file/format used by Pipe0_PreProcess
    struct ZDriftCut
    {
        double zMin;
        double zMax;
    };
    const double nSigmaZDrift {5.0};
    std::map<int, ZDriftCut> zDriftCuts;
    {
        std::ifstream finZ("../../Macros/L1PID/Outputs/zDrift_cut_perRun.dat");
        if(!finZ.is_open())
            throw std::runtime_error("Could not open zDrift_cut_perRun.dat");
        std::string line;
        while(std::getline(finZ, line))
        {
            if(line.empty() || line[0] == '#')
                continue;
            std::istringstream iss(line);
            int run;
            double meanS, sigmaS, meanE, sigmaE;
            if(iss >> run >> meanS >> sigmaS >> meanE >> sigmaE)
                zDriftCuts[run] = {meanS + nSigmaZDrift * sigmaS, meanE - nSigmaZDrift * sigmaE};
        }
    }

    // Select ammount of runs
    // dataManager.SetRuns(69, 70);

    // Get df for the runs
    auto chain {dataManager.GetChain(ActRoot::ModeType::EReadSilMod)};
    auto chain1 {dataManager.GetChain(ActRoot::ModeType::EMerge)};
    auto chain2 {dataManager.GetChain(ActRoot::ModeType::EFilter)};
    chain->AddFriend(chain1.get());
    chain->AddFriend(chain2.get(), "TPCData");
    auto df {ROOT::RDataFrame(*chain)};

    auto def {df.Filter([](ActRoot::ModularData& m, ActRoot::MergerData& mer) { return m.Get("GATCONF") == 4; },
                        {"ModularData", "MergerData"})}; // only f0

    // Print the number of entries in the dataframe
    auto nEntries = static_cast<int>(def.Count().GetValue());
    std::cout << "Number of entries in the dataframe: " << nEntries << std::endl;

    auto defFiltered = def.Filter(
        [validationZoneL1, driftFactor, zDriftCuts](ActRoot::TPCData& tpc, ActRoot::MergerData& m)
        {
            // First check for binary reactions
            auto sizeClusters = tpc.fClusters.size();
            if(sizeClusters != 3)
                return false;

            // Then check if light particle stops inside ACTAR

            // First check for the light cluster, which is the one with the largest angle with respect to the beam
            // direction
            int lightIdx {0};
            double maxAngle {-1111};
            ActRoot::Cluster lightCluster {};
            ROOT::Math::XYZVectorF beamlikeDir {1, 0, 0};
            for(int i = 0; i < sizeClusters; ++i)
            {
                const auto& cluster = tpc.fClusters[i];
                double angle = std::abs(TMath::ACos(cluster.GetLine().GetDirection().Unit().Dot(beamlikeDir)));
                if(angle > maxAngle)
                {
                    maxAngle = angle;
                    lightCluster = cluster;
                    lightIdx = i;
                }
            }
            // Now check if the cluster stops inside ACTAR
            auto lightLine = lightCluster.GetLine();
            lightCluster.SortAlongDir(lightLine.GetDirection());
            auto lastVoxel = lightCluster.GetVoxels().back();
            bool stopsInsideACTAR = (lastVoxel.GetPosition().X() < 128 - validationZoneL1 &&
                                     lastVoxel.GetPosition().X() > validationZoneL1 &&
                                     lastVoxel.GetPosition().Y() < 128 - validationZoneL1 &&
                                     lastVoxel.GetPosition().Y() > validationZoneL1);
            if(!stopsInsideACTAR)
                return false;

            // Now filter by charge and zDrift

            // start defining zDrift
            if(tpc.fRPs.empty())
                return false;
            auto rpVox = tpc.fRPs.front();
            auto& voxels = lightCluster.GetRefToVoxels();

            float maxDist = -1.0;
            float zExtreme = 0.0;
            for(auto& v : voxels)
            {
                auto pos = v.GetPosition();
                float dx = pos.X() - rpVox.X();
                float dy = pos.Y() - rpVox.Y();
                float dz = pos.Z() - rpVox.Z();
                float dist = std::sqrt(dx * dx + dy * dy + dz * dz);
                if(dist > maxDist)
                {
                    maxDist = dist;
                    zExtreme = pos.Z();
                }
            }
            float deltaZ = (zExtreme - rpVox.Z()) * driftFactor;
            auto zFinalPosition = 110.0f + deltaZ;
            // Now check if the zDrift is within the allowed range for the run
            auto it = zDriftCuts.find(m.fRun);
            if(it == zDriftCuts.end())
                return false; // run not in file → reject
            auto zDriftCut = it->second;
            if(zFinalPosition < zDriftCut.zMin || zFinalPosition > zDriftCut.zMax)
                return false;

            // Now the Qave
            auto cluster = lightCluster; // copy, as ComputeQave does
            cluster.ScaleVoxels(2, driftFactor);

            double gapThresh {6.5};
            double dist {};
            auto& voxelsForQave = cluster.GetVoxels();
            for(std::size_t v = 1; v < voxelsForQave.size(); v++)
            {
                auto p0 = cluster.GetLine().ProjectionPointOnLine(voxelsForQave[v - 1].GetPosition());
                auto p1 = cluster.GetLine().ProjectionPointOnLine(voxelsForQave[v].GetPosition());
                auto d = (p1 - p0).R();
                if(d < gapThresh)
                    dist += d;
            }
            auto qTotal = std::accumulate(voxelsForQave.begin(), voxelsForQave.end(), 0.f,
                                          [](float sum, const ActRoot::Voxel& v) { return sum + v.GetCharge(); });
            if(dist <= 0)
                return false;
            auto qAve = (qTotal / dist);
            if(qAve < 500)
                return false;

            return true;
        },
        {"TPCData", "MergerData"});


    // Get counts that stops inside ACTAR
    auto nEntriesFiltered = static_cast<int>(defFiltered.Count().GetValue());
    std::cout << "Number of entries that stops inside ACTAR: " << nEntriesFiltered << std::endl;

    // Save the filtered events to a root file
    std::cout << "Saving filtered events to file: ./Outputs/eventsL1_TriggerF0.dat" << std::endl;
    defFiltered.Snapshot("Filtered_Tree", "./Outputs/eventsL1_TriggerF0.root");

    // std::ofstream outFile("./Outputs/eventsL1_TriggerF0.dat");
    // defFiltered.Foreach([&outFile](ActRoot::MergerData& m) { m.Stream(outFile); }, {"MergerData"});
}