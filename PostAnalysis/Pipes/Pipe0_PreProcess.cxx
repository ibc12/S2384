#ifndef PIPE0_PREPROCESS_H
#define PIPE0_PREPROCESS_H
#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActParticle.h"
#include "ActSRIM.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TH1.h"
#include "TMath.h"
#include "TString.h"

#include <fstream>
#include <iostream>
#include <string>

#include "../../PrettyStyle.C"
#include "../HistConfig.h"

void Pipe0_PreProcess(const std::string& beam, const std::string& target, const std::string& light, bool isFiltered)
{
    ///////////////////////////////////////////////////////
    // Filter all the 7Li or the 11Li before doing the PID
    ///////////////////////////////////////////////////////

    std::string dataconf {};
    if(beam == "11Li")
        dataconf = "./../configs/data_11Li.conf";
    else if(beam == "7Li")
        dataconf = "./../configs/data_7Li.conf";
    else
        throw std::runtime_error("Beam cannot differ from 11Li or 7Li");

    // Read data
    ActRoot::DataManager dataman {dataconf, ActRoot::ModeType::EMerge};
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    auto chain3 {dataman.GetChain(ActRoot::ModeType::EFilter)};
    auto chain4 {dataman.GetChain(ActRoot::ModeType::EReadTPC)};
    chain->AddFriend(chain2.get());
    chain->AddFriend(chain3.get(), "TPCData");
    chain->AddFriend(chain4.get(), "GETTree");

    // Get drift parameter
    ActRoot::InputParser parser {};
    parser.ReadFile("../configs/detector.conf");
    auto driftBlock = parser.GetBlock("Merger");
    auto driftFactor = driftBlock->GetDouble("DriftFactor"); // in mm^2/us

    // Read per-run zDrift cuts: columns are run, meanStart, sigmaStart, meanEnd, sigmaEnd
    struct ZDriftCut
    {
        double zMin; // meanStart + nSigma * sigmaStart
        double zMax; // meanEnd   - nSigma * sigmaEnd
    };
    const double nSigmaZDrift {3.0};
    std::map<int, ZDriftCut> zDriftCuts;
    {
        std::ifstream finZ("../Macros/L1PID/Outputs/zDrift_cut_perRun.dat");
        if(!finZ.is_open())
            throw std::runtime_error("Could not open zDrift_cut_perRun.dat");
        std::string line;
        while(std::getline(finZ, line))
        {
            if(line.empty() || line[0] == '#')
                continue;
            std::istringstream iss(line);
            int run;
            double meanS, sigmaS, meanE, sigmaE;
            if(iss >> run >> meanS >> sigmaS >> meanE >> sigmaE)
                zDriftCuts[run] = {meanS + nSigmaZDrift * sigmaS, meanE - nSigmaZDrift * sigmaE};
        }
        finZ.close();
    }

    // RDataFrame
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {*chain};

    // 1-> Define zDrift for all events (maybe better only for L1)
    auto dfZ = df.Define("zDrift",
                         [&driftFactor](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                         {
                             auto& rp {m.fRP};

                             auto idx = m.fLightIdx;
                             if(idx < 0)
                                 return -111.0f;

                             if(static_cast<std::size_t>(idx) >= tpc.fClusters.size())
                                 return -111.0f;

                             auto& cluster = tpc.fClusters.at(idx);
                             auto& voxels = cluster.GetRefToVoxels();

                             if(tpc.fRPs.empty())
                                 return -111.0f;
                             auto rpVox = tpc.fRPs.front();

                             float maxDist = -1.0;
                             float zExtreme = 0.0;

                             for(auto& v : voxels)
                             {
                                 auto pos = v.GetPosition();

                                 float dx = pos.X() - rpVox.X();
                                 float dy = pos.Y() - rpVox.Y();
                                 float dz = pos.Z() - rpVox.Z();

                                 float dist = std::sqrt(dx * dx + dy * dy + dz * dz);

                                 if(dist > maxDist)
                                 {
                                     maxDist = dist;
                                     zExtreme = pos.Z();
                                 }
                             }

                             // diferencia en Z respecto al RP voxel
                             float deltaZ = (zExtreme - rpVox.Z()) * driftFactor;

                             float zDrift = 110.0 + deltaZ; // Force rp to be at 110 mm

                             return zDrift;
                         },
                         {"MergerData", "TPCData"});

    if(!isFiltered)
    {
        auto name {TString::Format("./Outputs/tree_preprocess_%s.root", beam.c_str())};
        std::cout << "Saving Preprocessed_Tree in file : " << name << '\n';
        dfZ.Snapshot("PreProcessed_Tree", name.Data());
        return;
    }

    // 1-> Filter the high charge deposits
    auto dfFilterCharge = dfZ.Filter( // Check for heavier clusters than Li
                                 [](ActRoot::MergerData& m)
                                 {
                                     if(m.fHeavy.fQave > 2000.)
                                         return false;
                                     return true;
                                 },
                                 {"MergerData"})
                              .Filter( // Check if heavy reaches end of ACTAR (other way to mask heavier clusters,
                                       // but can delete good events with heavy particle bad reconstructed)
                                  [](ActRoot::MergerData& m)
                                  {
                                      auto TL {m.fHeavy.fTL};
                                      auto theta {m.fThetaHeavy * TMath::DegToRad()};
                                      auto z_end {m.fRP.X() + TL * TMath::Cos(theta)};
                                      if(z_end < 240)
                                          return false;
                                      return true;
                                  },
                                  {"MergerData"})
                              .Filter( // Most of the times high charge deposit
                                       // is masked by rp or is in beam cluster
                                  [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                                  {
                                      auto rp {m.fRP};
                                      auto rp_y {rp.Y() / 2};
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
                                  {"MergerData", "GETTree.TPCData"});


    // 3-> Filter events with bad zDrift and low charge (only for L1 events)
    auto dfFilterZandQ = dfFilterCharge.Filter(
        [&zDriftCuts](float zDrift, ActRoot::MergerData& m, ActRoot::ModularData& mod)
        {
            if(mod.Get("GATCONF") != 8)
                return true; // not L1 -> no problem going to pad plane / cathode
            auto it = zDriftCuts.find(m.fRun);
            if(it == zDriftCuts.end())
                return false; // run not in file → reject
            return zDrift > it->second.zMin && zDrift < it->second.zMax && m.fLight.fQave > 600;
        },
        {"zDrift", "MergerData", "ModularData"});

    if(isFiltered)
    {
        auto name {TString::Format("./Outputs/tree_preprocess_F_%s.root", beam.c_str())};
        std::cout << "Saving PreProcessed_Tree in file : " << name << '\n';
        dfFilterZandQ.Snapshot("PreProcessed_Tree", name.Data());
    }
}

#endif