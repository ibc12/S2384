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
    for(const auto& [layer, n] : it->second)
    {
        auto itLayer = std::find(layers.begin(), layers.end(), layer);
        if(itLayer != layers.end())
        {
            auto idx = std::distance(layers.begin(), itLayer);
            auto itVec = fSilData->fSiN.find(layer);
            if(itVec == fSilData->fSiN.end())
                continue;
            const auto& vec = itVec->second;
            if(idx < vec.size() && vec[idx] == n)
            {
                // Neutralize that silicon hit by setting its calibrated energy to 0
                // later threshold checks will ignore the hit because energy is zero.
                if(idx < fSilData->fSiE[layer].size())
                {
                    fSilData->fSiE[layer][idx] = 0.0f;
                }
            }
        }
    }
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