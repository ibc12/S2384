#include "ActVAction.h"

#include "ActTPCParameters.h"
#include "ActContinuity.h"

namespace ActAlgorithm
{
class BrokenTracksZ : public VAction
{
public:
    // Parameters of the action
    int fRebinX {};
    int fRebinY {};
    int fRebinZ {};
    int fVoxelsContinuity {};

    ActRoot::TPCParameters fTPCParsRebined {};
    std::shared_ptr<ActAlgorithm::Continuity> fContinuity {};
public:
    BrokenTracksZ() : VAction("BrokenTracksZ") {}

    void ReadConfiguration(std::shared_ptr<ActRoot::InputBlock> block) override;
    void Run() override;
    void Print() const override;
};
} // namespace ActAlgorithm