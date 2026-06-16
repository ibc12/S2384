#ifndef Task2_h
#define Task2_h

#include "ActVTask.h"

#include <map>
#include <string>
#include <vector>

// Task2: Erase multiplicity in some silicon layers if there is any.
// For example F2 events with multiplicity > 1 will be erased. Posible noise, also  to avoid pileup events (not able to
// distinguish which is the real hit).

namespace ActAlgorithm
{
class Task2 : public VTask
{
public:
    std::vector<std::string> fLayerToErase {"f2"}; // Layer to erase multiplicity

    Task2() : VTask("Task2") {}

    bool Run() override;
    void Print() override;
};
} // namespace ActAlgorithm

#endif // !Task2_h