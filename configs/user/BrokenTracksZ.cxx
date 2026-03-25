#include "BrokenTracksZ.h"

#include "ActAFindRP.h"
#include "ActAlgoFuncs.h"
#include "ActCluster.h"
#include "ActColors.h"
#include "ActContinuity.h"
#include "ActInputParser.h"
#include "ActMergerData.h"
#include "ActMultiAction.h"
#include "ActRANSAC.h"
#include "ActTPCData.h"

#include "TMath.h"

#include <algorithm>
#include <memory>

void ActAlgorithm::BrokenTracksZ::ReadConfiguration(std::shared_ptr<ActRoot::InputBlock> block)
{
    fIsEnabled = block->GetBool("IsEnabled");
    // if(!fIsEnabled) // Comment if putting in FinerActions Isenabled forced true after parsing
    //     return;
    if(block->CheckTokenExists("rebinX"))
        fRebinX = block->GetInt("rebinX");
    if(block->CheckTokenExists("rebinY"))
        fRebinY = block->GetInt("rebinY");
    if(block->CheckTokenExists("rebinZ"))
        fRebinZ = block->GetInt("rebinZ");
    if(block->CheckTokenExists("voxelsContinuity"))
        fVoxelsContinuity = block->GetInt("voxelsContinuity");
}

void ActAlgorithm::BrokenTracksZ::Run()
{
    if(!fIsEnabled)
        return;
    // We cannot set tpcPars and continuity in Readconfiguration, because fTPCPars is not seted yet
    if(!fContinuity)
    {
        fTPCParsRebined = ActRoot::TPCParameters {
            (fTPCPars->GetNPADSX() / fRebinX),
            (fTPCPars->GetNPADSY() / fRebinY),
            (fTPCPars->GetNPADSZ() / fRebinZ),
        };
        fContinuity = std::make_shared<ActAlgorithm::Continuity>(&fTPCParsRebined, fVoxelsContinuity);
    }
    // Get noise and light cluster
    const auto& noise {fTPCData->fRaw};
    // Get particle with higher angle respect (1,0,0) -> LightCluster and save position
    int lightIdx {0};
    double maxAngle {-1111};
    ActRoot::Cluster lightCluster {};
    ROOT::Math::XYZVectorF beamlikeDir {1, 0, 0};
    for(int i = 0; i < fTPCData->fClusters.size(); ++i)
    {
        const auto& cluster = fTPCData->fClusters[i];
        double angle = std::abs(TMath::ACos(cluster.GetLine().GetDirection().Unit().Dot(beamlikeDir)));
        if(angle > maxAngle)
        {
            maxAngle = angle;
            lightCluster = cluster;
            lightIdx = i;
        }
    }
    // Join all voxels of noise and light in a single vector to do the rebinning and continuity
    std::vector<ActRoot::Voxel> voxels;
    voxels.insert(voxels.end(), lightCluster.GetVoxels().begin(), lightCluster.GetVoxels().end());
    voxels.insert(voxels.end(), noise.begin(), noise.end());

    // Run the rebining
    auto result = RebinTracks(voxels, fTPCPars, fRebinX, fRebinY, fRebinZ);

    auto clusterRebinedAfterContinuity = fContinuity->Run(result.first, true);
    // Get bigger cluster after continuity
    std::vector<ActRoot::Voxel> VoxelsToUnrebin;
    if(!clusterRebinedAfterContinuity.first.empty())
    {
        auto maxCluster =
            std::max_element(clusterRebinedAfterContinuity.first.begin(), clusterRebinedAfterContinuity.first.end(),
                             [](const ActRoot::Cluster& a, const ActRoot::Cluster& b)
                             { return a.GetVoxels().size() < b.GetVoxels().size(); });
        VoxelsToUnrebin = maxCluster->GetRefToVoxels();
    }
    else
    {
        if(fIsVerbose)
        {
            std::cout << "No clusters found after continuity, skipping event" << '\n';
        }
        return;
    }

    // Undo rebinning
    auto unrebinedClusterAfterContinuity = UndoRebinning(VoxelsToUnrebin, result.second);
    //
    if(unrebinedClusterAfterContinuity.size() < lightCluster.GetVoxels().size())
    {
        if(fIsVerbose)
        {
            std::cout << BOLDRED << "Continuity broke the track, skipping event" << RESET << '\n';
        }
        return;
    }
    if(unrebinedClusterAfterContinuity.size() > lightCluster.GetVoxels().size())
    {
        if(fIsVerbose)
        {
            std::cout << BOLDGREEN << "Continuity added voxels to the track, do substitution of cluster" << RESET
                      << '\n';
            std::cout << BOLDGREEN << "Light cluster before processing: " << lightCluster.GetVoxels().size()
                      << " voxels" << '\n';
            std::cout << BOLDGREEN << "Cluster after processing: " << unrebinedClusterAfterContinuity.size()
                      << " voxels" << '\n';

            std::cout << BOLDGREEN << "Diferent voxels positions: " << '\n';
            // Go through the after voxels, and compare the positions with all the lightClustre voxels
            for(auto& voxelAfter : unrebinedClusterAfterContinuity)
            {
                auto posAfter = voxelAfter.GetPosition();
                bool found = false;
                for(auto& voxelBefore : lightCluster.GetVoxels())
                {
                    auto posBefore = voxelBefore.GetPosition();
                    if(posAfter == posBefore)
                    {
                        found = true;
                        break;
                    }
                }
                if(!found)
                {
                    std::cout << "Voxel added by continuity: " << posAfter.X() << " " << posAfter.Y() << " "
                              << posAfter.Z() << '\n';
                }
            }
        }
        // Erase the lightCluster index form fTPCData->fClusters and add the unrebinedClusterAfterContinuity as a new
        // cluster
        ActRoot::Cluster newCluster {};
        newCluster.SetVoxels(unrebinedClusterAfterContinuity);
        newCluster.ReFit();
        newCluster.ReFillSets();
        auto& line = newCluster.GetLine(); // Force to not have nans in the directions
        if(std::isnan(line.GetDirection().X()) || std::isnan(line.GetDirection().Y()) || std::isnan(line.GetDirection().Z()))
        {
            if(fIsVerbose)
            {
                std::cout << BOLDRED << "Continuity added voxels but the fit failed, skipping event" << RESET << '\n';
            }
            return;
        }
        fTPCData->fClusters[lightIdx] = newCluster;
        return;
    }
    if(unrebinedClusterAfterContinuity.size() == lightCluster.GetVoxels().size())
    {
        if(fIsVerbose)
        {
            std::cout << BOLDBLUE << "Continuity kept the track, skipping event" << RESET << '\n';
        }
        return;
    }
}

void ActAlgorithm::BrokenTracksZ::Print() const
{
    std::cout << BOLDCYAN << "····· " << GetActionID() << " ·····" << '\n';
    if(!fIsEnabled)
    {
        std::cout << "····························" << RESET << '\n';
        return;
    }
    std::cout << "RebinX:           " << fRebinX << RESET << '\n';
    std::cout << "RebinY:           " << fRebinY << RESET << '\n';
    std::cout << "RebinZ:           " << fRebinZ << RESET << '\n';
    std::cout << "VoxelsContinuity: " << fVoxelsContinuity << RESET << '\n';
}

// Create symbol to load class from .so
extern "C" ActAlgorithm::BrokenTracksZ* CreateUserAction()
{
    return new ActAlgorithm::BrokenTracksZ;
}
