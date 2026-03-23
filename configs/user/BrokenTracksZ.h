#include "ActVAction.h"

namespace ActAlgorithm
{
class BrokenTracksZ : public VAction
{
public:
    // Parameters of the action
public:
    BrokenTracksZ() : VAction("BrokenTracksZ") {}

    void ReadConfiguration(std::shared_ptr<ActRoot::InputBlock> block) override;
    void Run() override;
    void Print() const override;
};
} // namespace ActAlgorithm