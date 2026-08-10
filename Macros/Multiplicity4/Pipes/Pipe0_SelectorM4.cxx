#ifndef Pipe0_SelectorM4_cxx
#define Pipe0_SelectorM4_cxx

#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActParticle.h"
#include "ActSRIM.h"
#include "ActSilData.h"
#include "ActSilMatrix.h"
#include "ActSilSpecs.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TH2.h"
#include "TH2D.h"
#include "TString.h"

#include <map>
#include <string>

#include "../../../PrettyStyle.C"
#include "../Utils.h"

void Pipe0_SelectorM4(const std::string& beam, const std::string& target, const std::string& light)
{
    // As I did it before I depend on the Merger to get the data, with that I cannot recover the f0 information, because
    // the merger erase events with multiplicity two in the slicons

    std::string dataconf {};
    if(beam == "11Li")
        dataconf = "./../../configs/data_11Li.conf";
    else if(beam == "7Li")
        dataconf = "./../../configs/data_7Li.conf";
    else
        throw std::runtime_error("Beam cannot differ from 11Li or 7Li");

    // Read data
    ActRoot::DataManager dataman {dataconf, ActRoot::ModeType::EFilter};
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    auto chain4 {dataman.GetChain(ActRoot::ModeType::EReadTPC)};
    auto chain5 {dataman.GetChain(ActRoot::ModeType::EMerge)};
    chain->AddFriend(chain2.get());
    chain->AddFriend(chain4.get(), "GETTree");
    chain->AddFriend(chain5.get());

    // Get drift parameter
    ActRoot::InputParser parser {};
    parser.ReadFile("../../configs/detector.conf");
    auto driftBlock = parser.GetBlock("Merger");
    auto driftFactor = driftBlock->GetDouble("DriftFactor"); // in mm^2/us

    // RDataFrame
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {*chain};

    // df.Describe().Print();

    // First filter TPC multiplicity 4 and ensure a later silicon hit
    auto dfFilter = df.Filter("fClusters.size() == 4")
                        .Filter(
                            [](ActRoot::ModularData& m)
                            {
                                return (m.Get("GATCONF") == 1 || m.Get("GATCONF") == 2 || m.Get("GATCONF") == 4 ||
                                        m.Get("GATCONF") == 8);
                            },
                            {"ModularData"}) // All GATCONF that can have events
                        .Filter(
                            [](ActRoot::SilData& sil, ActRoot::ModularData& m)
                            {
                                if(m.Get("GATCONF") == 1)
                                {
                                    if(sil.fSiN.at("l0").front() == 9) // Filter l0_9, there's no matrix for that
                                        return false;
                                    else
                                        return true;
                                }
                                if(m.Get("GATCONF") == 4)
                                {
                                    if(sil.fSiN.at("f0").front() == 10 || sil.fSiN.at("f0").front() == 9 ||
                                       sil.fSiN.at("f0").front() == 5) // there's no matrix for that
                                        return false;
                                    else
                                        return true;
                                }
                                else
                                    return true;
                            },
                            {"SilData", "ModularData"})
                        .Filter( // Most of the times high charge deposit
                                 // is masked by rp or is in beam cluster
                            [&](ActRoot::TPCData& f, ActRoot::TPCData& tpc)
                            {
                                // Check first if the RP is defined, otherwise return false
                                if(f.fRPs.empty())
                                    return false;
                                auto rp {f.fRPs.front()};
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
                                           v.GetPosition().X() > rp_x - 10 && v.GetPosition().X() < rp_x + 10)
                                            if(v.GetCharge() > 3000.)
                                                counter++;
                                    }
                                }
                                if(counter > 4)
                                    return false;
                                return true;
                            },
                            {"TPCData", "GETTree.TPCData"})
                        .Filter(
                            [](ActRoot::TPCData& tpc)
                            {
                                // Ensure there is one beam-like cluster
                                int bl_counter = 0;
                                for(auto& cluster : tpc.fClusters)
                                {
                                    if(cluster.GetIsBeamLike())
                                        bl_counter++;
                                }
                                return bl_counter == 1;
                            },
                            {"TPCData"});

    auto dfLight = dfFilter.Define("LightIdx", [&](ActRoot::TPCData& tpc)
                                       { // Get index of light particle - highest angle with respect to beam-like particle
        auto beamIdx = -1;
        for(int i = 0; i < tpc.fClusters.size(); ++i)
        {
            auto& cluster = tpc.fClusters[i];
            if(cluster.GetIsBeamLike())
            {
                beamIdx = i;
                break;;
            }
        }
        if(beamIdx == -1)
            return -1;
        auto& beamCluster = tpc.fClusters[beamIdx];
        auto beamDir = beamCluster.GetLine().GetDirection();
        double maxAngle = -1.0;
        int lightIdx = -1;
        for(int i = 0; i < tpc.fClusters.size(); ++i)
        {
            if(i == beamIdx)
                continue;
            auto& cluster = tpc.fClusters[i];
            auto clusterDir = cluster.GetLine().GetDirection();
            auto angle = Utils::GetTheta3D(beamDir, clusterDir);
            if(angle > maxAngle)
            {
                maxAngle = angle;
                lightIdx = i;
            }
        }
        return lightIdx;
        }, {"TPCData"})
            .Define("zDrift",
            [&driftFactor](ActRoot::TPCData& tpc, int idx)
            {
                auto& cluster = tpc.fClusters.at(idx);
                auto& voxels = cluster.GetRefToVoxels();

                if(tpc.fRPs.empty())
                    return -111.0f;
                auto rpVox = tpc.fRPs.front();

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
                // diferencia en Z respecto al RP voxel
                float deltaZ = (zExtreme - rpVox.Z()) * driftFactor;

                float zDrift = 110.0 + deltaZ; // Force rp to be at 110 mm

                return zDrift;
            },
            {"TPCData", "LightIdx"})
        .Define("Qave",
                [](ActRoot::TPCData& tpc, int idx)
                {
                    auto cluster = tpc.fClusters[idx];
                    cluster.SortAlongDir(cluster.GetLine().GetDirection());
                    auto voxels = cluster.GetRefToVoxels();
                    float Q = 0.f;
                    for(const auto& v : voxels)
                        Q += v.GetCharge();
                    ROOT::Math::XYZPointF firstPos = voxels.front().GetPosition();
                    ROOT::Math::XYZPointF lastPos = voxels.back().GetPosition();
                    Utils::ScalePoint(firstPos);
                    Utils::ScalePoint(lastPos);
                    float length = (lastPos - firstPos).R();
                    if(length <= 0.0)
                        return -1.f;
                    return Q / length;
                },
                {"TPCData", "LightIdx"});

    std::cout << "Counts after filtering: " << dfFilter.Count().GetValue() << std::endl;

    // Save dataframe in a .root file
    TString outfile = TString::Format("./Outputs/SelectorM4_%s.root", beam.c_str());
    dfLight.Snapshot("Final_Tree", outfile);
    std::cout << "Saving Final_Tree in " << outfile << " from Pipe0_SelectorM4" << '\n';
}
#endif