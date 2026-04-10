#include "ActVAction.h"

namespace ActAlgorithm
{
class PrintAfterFindRP : public VAction
{
public:
    // Parameters of the action
public:
    PrintAfterFindRP() : VAction("PrintAfterFindRP") {}

    void ReadConfiguration(std::shared_ptr<ActRoot::InputBlock> block) override;
    void Run() override;
    void Print() const override;
};
} // namespace ActAlgorithm
