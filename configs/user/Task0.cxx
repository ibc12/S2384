#include "Task0.h"

#include "ActMergerData.h"
#include "ActSilData.h"
#include "ActVTask.h"

#include <iostream>


bool ActAlgorithm::Task0::Run()
{
    // If no silicon impact, do nothing
    if(fSilData->GetLayers().empty())
        return true;

    auto it = fMap.find(fMergerData->fRun);
    if(it == fMap.end())
        return true; // Run not affected by any silicon issue, do nothing

    const auto& layers = fSilData->GetLayers();
    for(const auto& [name, n] : it->second)
    {
        auto itLayer = std::find(layers.begin(), layers.end(), name);
        if(itLayer != layers.end())
        {
            // Get layer
            auto& layerNs = fSilData->fSiN[name];
            auto itN = std::find(layerNs.begin(), layerNs.end(), n);
            if(itN == layerNs.end())
                continue; // This silicon hit is not present in this event, do nothing
            auto idx = std::distance(layerNs.begin(), itN);
            // Neutralize that silicon hit by setting its calibrated energy to 0
            // later threshold checks will ignore the hit because energy is zero.
            if(idx < fSilData->fSiE[name].size())
            {
                fSilData->fSiE[name][idx] = 0.0f;
            }
        }
    }
    // fSilData->Print(); // Print the silicon data for this event to check that the changes were applied correctly
    return true;
}

void ActAlgorithm::Task0::Print()
{
    std::cout << "Task0 prints!" << '\n';
}

// Create symbol to load class from .so
// extern "C" disables C++ function name mangling
extern "C" ActAlgorithm::Task0* Create()
{
    return new ActAlgorithm::Task0;
}