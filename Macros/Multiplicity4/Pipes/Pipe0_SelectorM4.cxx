#ifndef Pipe0_SelectorM4_cxx
#define Pipe0_SelectorM4_cxx

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

#include <map>
#include <string>

#include "../../../PrettyStyle.C"
#include "../Utils.h"

void Pipe0_SelectorM4(const std::string& beam, const std::string& target, const std::string& light)
{
    // As I did it before I depend on the Merger to get the data, with that I cannot recover the f0 information, because
    // the merger erase events with multiplicity two in the slicons

    std::string dataconf {};
    if(beam == "11Li")
        dataconf = "./../../configs/data_11Li.conf";
    else if(beam == "7Li")
        dataconf = "./../../configs/data_7Li.conf";
    else
        throw std::runtime_error("Beam cannot differ from 11Li or 7Li");

    // Read data
    ActRoot::DataManager dataman {dataconf, ActRoot::ModeType::EFilter};
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    auto chain4 {dataman.GetChain(ActRoot::ModeType::EReadTPC)};
    auto chain5 {dataman.GetChain(ActRoot::ModeType::EMerge)};
    chain->AddFriend(chain2.get());
    chain->AddFriend(chain4.get(), "GETTree");
    chain->AddFriend(chain5.get());

    // RDataFrame
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {*chain};

    // df.Describe().Print();

    // First filter TPC multiplicity 4 and ensure a later silicon hit
    auto dfFilter = df.Filter("fClusters.size() == 4")
                        .Filter([](ActRoot::ModularData& m)
                                { return (m.Get("GATCONF") == 1 || m.Get("GATCONF") == 2); }, {"ModularData"})
                        .Filter(
                            [](ActRoot::SilData& sil, ActRoot::ModularData& m)
                            {
                                if(m.Get("GATCONF") == 1)
                                    if(sil.fSiN.at("l0").front() == 9) // Filter l0_9, there's no matrix for that
                                        return false;
                                    else
                                        return true;
                                else
                                    return true;
                            },
                            {"SilData", "ModularData"})
                        .Filter( // Most of the times high charge deposit
                                 // is masked by rp or is in beam cluster
                            [&](ActRoot::TPCData& f, ActRoot::TPCData& tpc)
                            {
                                auto rp {f.fRPs.front()};
                                auto rp_y {rp.Y()};
                                // Run for all clusters
                                int counter {};
                                for(auto& cluster : tpc.fClusters)
                                {
                                    auto voxels {cluster.GetRefToVoxels()};
                                    for(auto& v : voxels)
                                    {
                                        if(v.GetPosition().Y() > rp_y - 3 &&
                                           v.GetPosition().Y() < rp_y + 3) // aprox L1 exclusion zone
                                            if(v.GetCharge() > 3000.)
                                                counter++;
                                    }
                                }
                                if(counter > 4)
                                    return false;
                                return true;
                            },
                            {"TPCData", "GETTree.TPCData"})
                        .Filter(
                            [](ActRoot::TPCData& tpc)
                            {
                                // Ensure there is one beam-like cluster
                                int bl_counter = 0;
                                for(auto& cluster : tpc.fClusters)
                                {
                                   if(cluster.GetIsBeamLike())
                                       bl_counter++;
                                }
                                return bl_counter == 1;
                            },
                            {"TPCData"});

    // Save dataframe in a .root file
    TString outfile = TString::Format("./Outputs/SelectorM4_%s.root", beam.c_str());
    dfFilter.Snapshot("Final_Tree", outfile);
    std::cout << "Saving Final_Tree in " << outfile << " from Pipe0_SelectorM4" << '\n';
}
#endif