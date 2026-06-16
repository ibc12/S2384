#include "Task2.h"

#include "ActMergerData.h"
#include "ActSilData.h"
#include "ActVTask.h"

#include <iostream>


bool ActAlgorithm::Task2::Run()
{
    // If no silicon impact, do nothing
    if(fSilData->GetLayers().empty())
        return true;
    for(auto& layer : fLayerToErase)
    {
        const auto& layers = fSilData->GetLayers();
        auto itLayer = std::find(layers.begin(), layers.end(), layer);
        if(itLayer != layers.end())
        {
            // Get layer
            auto& layerNs = fSilData->fSiN[layer];
            if(layerNs.size() > 1)
            {
                // Erase multiplicity by keeping only the first hit and neutralizing the rest
                for(int i = 0; i < layerNs.size(); ++i)
                {
                    fSilData->fSiE[layer][i] = 0.0f; // Neutralize the hit by setting its calibrated energy to 0
                }
            }
        }
    }
    return true;
}

void ActAlgorithm::Task2::Print()
{
    std::cout << "Task2 prints!" << '\n';
    std::cout << "Layers to erase multiplicity: ";
    for(auto& layer : fLayerToErase)
        std::cout << layer << " ";
    std::cout << '\n';
}

// Create symbol to load class from .so
// extern "C" disables C++ function name mangling
extern "C" ActAlgorithm::Task2* Create()
{
    return new ActAlgorithm::Task2;
}