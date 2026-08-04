#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActParticle.h"
#include "ActSRIM.h"
#include "ActSilData.h"
#include "ActSilMatrix.h"
#include "ActSilSpecs.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TH2.h"
#include "TH2D.h"
#include "TString.h"

#include <iostream>
#include <map>
#include <string>

#include "../../../PrettyStyle.C"
#include "../Utils.h"


void checkL1NotValEventWithSilInfo()
{
    // Get Output from pipe0
    ROOT::RDataFrame df("Final_Tree", TString::Format("../Outputs/SelectorM4_7Li.root").Data());

    auto sils {std::make_shared<ActPhysics::SilSpecs>()};
    sils->ReadFile("../../../configs/silspecs.conf");

    // Event r72e85176 seems to have silicon information although the merger output is L1 not val so let's get the
    // information in this event

    df.Filter([](ActRoot::MergerData m) { return m.fRun == 72 && m.fEntry == 85176; }, {"MergerData"})
        .Foreach(
            [&](ActRoot::MergerData m, ActRoot::TPCData tpc, int LightIdx, ActRoot::SilData silIn,
                ActRoot::ModularData modIn)
            {
                auto gat = modIn.Get("GATCONF");
                std::cout << "GATCONF = " << gat << std::endl;

                ActRoot::SilData sil {silIn};
                sil.ApplyFinerThresholds(sils);

                // Get number of cluster
                std::cout << "Number of clusters in TPC: " << tpc.fClusters.size() << std::endl;

                // Get hist in each layer
                for(const std::string layer : {"l0", "r0", "f0"})
                {
                    auto it = sil.fSiN.find(layer);
                    if(it == sil.fSiN.end() || it->second.empty())
                    {
                        std::cout << "No hit in layer " << layer << std::endl;
                        continue;
                    }
                    auto hits = it->second;
                    std::cout << "Hits in layer " << layer << ": ";
                    for(int i = 0; i < hits.size(); ++i)
                    {
                        std::cout << "Idx hit: " << hits[i] << " Energy: " << sil.fSiE.at(layer)[i] << std::endl;
                    }
                }

                m.Print();
            },
            {"MergerData", "TPCData", "LightIdx", "SilData", "ModularData"});
}