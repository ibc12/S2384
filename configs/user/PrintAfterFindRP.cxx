#include "PrintAfterFindRP.h"

#include "ActAFindRP.h"
#include "ActCluster.h"
#include "ActColors.h"
#include "ActInputParser.h"
#include "ActMultiAction.h"
#include "ActRANSAC.h"
#include "ActTPCData.h"

#include <algorithm>
#include <memory>

void ActAlgorithm::PrintAfterFindRP::ReadConfiguration(std::shared_ptr<ActRoot::InputBlock> block)
{
    fIsEnabled = block->GetBool("IsEnabled");
    if(!fIsEnabled)
        return;
    // if(block->CheckTokenExists("MaxAngle"))
    //     fMaxAngle = block->GetDouble("MaxAngle");
    // if(block->CheckTokenExists("MinLength"))
    //     fMinLength = block->GetDouble("MinLength");
    // if(block->CheckTokenExists("CylinderR"))
    //     fCylinderR = block->GetDouble("CylinderR");
}

//////////////////////////////////////
// Action to print voxel information after FindRP, to check if voxels are sorted
//////////////////////////////////////
void ActAlgorithm::PrintAfterFindRP::Run()
{
    if(!fIsEnabled)
        return;

    if(fIsVerbose)
    {
        for(auto cluster : fTPCData->fClusters)
        {
            auto line = cluster.GetLine();
            cluster.SortAlongDir(line.GetDirection().Unit());

            std::cout << "======================" << '\n';
            std::cout << "Cluster with " << cluster.GetVoxels().size() << " voxels" << '\n';
            for(auto voxel : cluster.GetVoxels())
            {
                auto pos = voxel.GetPosition();
                std::cout << "Voxel position: " << pos.X() << " " << pos.Y() << " " << pos.Z() << '\n';
            }
        }
    }
}

void ActAlgorithm::PrintAfterFindRP::Print() const
{
    std::cout << BOLDCYAN << "····· " << GetActionID() << " ·····" << '\n';
    if(!fIsEnabled)
    {
        std::cout << "······························" << RESET << '\n';
        return;
    }
    std::cout << "  No parameters yet" << RESET << '\n';
}

// Create symbol to load class from .so
extern "C" ActAlgorithm::PrintAfterFindRP* CreateUserAction()
{
    return new ActAlgorithm::PrintAfterFindRP;
}
