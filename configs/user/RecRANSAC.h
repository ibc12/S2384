#include "ActVAction.h"

namespace ActAlgorithm
{
class RecRANSAC : public VAction
{
public:
    // Parameters of the action
public:
    RecRANSAC() : VAction("RecRANSAC") {}

    void ReadConfiguration(std::shared_ptr<ActRoot::InputBlock> block) override;
    void Run() override;
    void Print() const override;
};
} // namespace ActAlgorithm
