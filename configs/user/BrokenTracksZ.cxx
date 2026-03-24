#include "BrokenTracksZ.h"

#include "ActAFindRP.h"
#include "ActCluster.h"
#include "ActColors.h"
#include "ActInputParser.h"
#include "ActMultiAction.h"
#include "ActRANSAC.h"
#include "ActTPCData.h"
#include "ActAlgoFuncs.h"

#include <algorithm>
#include <memory>

void ActAlgorithm::BrokenTracksZ::ReadConfiguration(std::shared_ptr<ActRoot::InputBlock> block)
{
    fIsEnabled = block->GetBool("IsEnabled");
    if(!fIsEnabled)
        return;
    // if(block->CheckTokenExists("MaxAngle"))
    //     fMaxAngle = block->GetDouble("MaxAngle");
}

void ActAlgorithm::BrokenTracksZ::Run()
{
    if(!fIsEnabled)
        return;

    // Get noise and light cluster
    const auto& noise {fTPCData->fRaw};

    
    if(fIsVerbose)
    {
        std::cout << BOLDGREEN << "-- BrokenTracksZ --" << '\n';
        std::cout << "Noise voxels before: " << noise.size() << '\n';
        std::cout << "Noise voxels appended: " << noise.size() << RESET << '\n';
    }
}

void ActAlgorithm::BrokenTracksZ::Print() const
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
extern "C" ActAlgorithm::BrokenTracksZ* CreateUserAction()
{
    return new ActAlgorithm::BrokenTracksZ;
}
