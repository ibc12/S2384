#include "RecRANSAC.h"

#include "ActAFindRP.h"
#include "ActCluster.h"
#include "ActColors.h"
#include "ActInputParser.h"
#include "ActMultiAction.h"
#include "ActRANSAC.h"
#include "ActTPCData.h"

#include <algorithm>
#include <memory>

void ActAlgorithm::RecRANSAC::ReadConfiguration(std::shared_ptr<ActRoot::InputBlock> block)
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

void ActAlgorithm::RecRANSAC::Run()
{
    if(!fIsEnabled)
        return;

    // Determine beam likes in event, calling DetermineBeamLikes from FindRP before
    if(fMultiAction->HasAction("FindRP"))
    {
        auto findrp {std::dynamic_pointer_cast<ActAlgorithm::Actions::FindRP>(fMultiAction->GetAction("FindRP"))};
        if(findrp)
            findrp->ExecInnerAction("DetermineBeamLikes");
    }

    // Trigger only when all are beam-likes
    auto trigger {std::all_of(fTPCData->fClusters.begin(), fTPCData->fClusters.end(),
                              [](const ActRoot::Cluster& cl) { return cl.GetIsBeamLike(); })};
    if(!trigger)
        return;
    // if(fTPCData->fClusters.size() > 1)
    //     return;

    const auto& noise {fTPCData->fRaw};
    int iter {175};
    int minVoxels {7};
    double distThresh {2.};
    ActAlgorithm::RANSAC ransac {iter, minVoxels, distThresh};
    auto [clusters, back] {ransac.Run(noise)};
    if(fIsVerbose)
    {
        std::cout << BOLDGREEN << "-- RecRANSAC --" << '\n';
        std::cout << "Noise size: " << noise.size() << '\n';
        std::cout << "New clusters: " << clusters.size() << RESET << '\n';
    }
    if(clusters.size())
    {
        // Get clusters with best chi2
        std::sort(clusters.begin(), clusters.end(), [](const ActRoot::Cluster& a, const ActRoot::Cluster& b)
                  { return a.GetLine().GetChi2() < b.GetLine().GetChi2(); });
        // Push the first one
        auto& cluster {clusters.front()};
        // Set flag
        cluster.SetFlag("IsRANSAC", true);
        // Push to vector
        fTPCData->fClusters.push_back(cluster);
    }
}

void ActAlgorithm::RecRANSAC::Print() const
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
extern "C" ActAlgorithm::RecRANSAC* CreateUserAction()
{
    return new ActAlgorithm::RecRANSAC;
}
