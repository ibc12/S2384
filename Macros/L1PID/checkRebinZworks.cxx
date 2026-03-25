#include "ActCluster.h"
#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActLine.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActSRIM.h"
#include "ActTPCData.h"
#include "ActTPCParameters.h"
#include "ActVoxel.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"
#include <random>

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMath.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TSpline.h>
#include <TStyle.h>

#include <Math/Point3D.h>
#include <Math/Vector3D.h>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

#include "../../PrettyStyle.C"


void checkRebinZworks()
{
    PrettyStyle(true, true);
    bool isDiffusion = true; // Change between diffusion spline or not

    // Get all L1 experimental data for 7Li; cut in L1 and filter bad events, select elastic, check good fits to
    // deuterium

    // Get all data for 7Li (there is no triton there, so easier to see deuterium)
    std::string beam {"7Li"};
    std::string target {"d"};
    std::string light {"p"};

    std::string dataconf {};
    std::string dataconfAux {"../../configs/data_7Li_auxiliar.conf"};
    if(beam == "11Li")
        dataconf = "../../configs/data_11Li.conf";
    else if(beam == "7Li")
        dataconf = "../../configs/data_7Li.conf";
    else
        throw std::runtime_error("Beam cannot differ from 11Li or 7Li");

    // Read data
    ActRoot::DataManager dataman {dataconf, ActRoot::ModeType::EMerge};
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    auto chain3 {dataman.GetChain(ActRoot::ModeType::EFilter)};
    chain->AddFriend(chain2.get());
    chain->AddFriend(chain3.get());

    ActRoot::DataManager datamanAux {dataconfAux, ActRoot::ModeType::EMerge};
    auto chainAux {datamanAux.GetChain()};
    auto chainAux2 {datamanAux.GetChain(ActRoot::ModeType::EReadSilMod)};
    auto chainAux3 {datamanAux.GetChain(ActRoot::ModeType::EFilter)};
    chainAux->AddFriend(chainAux2.get());
    chainAux->AddFriend(chainAux3.get());

    // Get drift parameter
    ActRoot::InputParser parser {};
    parser.ReadFile("../../configs/detector.conf");
    auto driftBlock = parser.GetBlock("Merger");
    auto driftFactor = driftBlock->GetDouble("DriftFactor"); // in mm^2/us

    // RDataFrame
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame dforigin {*chain};
    ROOT::RDataFrame dforiginAux {*chainAux};

    // Filter GATCONF == L1
    auto df =
        dforigin.Filter([](ActRoot::ModularData& mod, ActRoot::MergerData& mer)
                        { return mod.Get("GATCONF") == 8 && (mer.fLightIdx != -1); }, {"ModularData", "MergerData"});

    auto dfAux =
        dforiginAux.Filter([](ActRoot::ModularData& mod, ActRoot::MergerData& mer)
                           { return mod.Get("GATCONF") == 8 && (mer.fLightIdx != -1); }, {"ModularData", "MergerData"});


    // Filter bad events
    auto dfFilter = df.Filter( // Check for heavier clusters than Li
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
                            [&](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
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
                            {"MergerData", "TPCData"});
    auto dfFilterAux = dfAux
                           .Filter( // Check for heavier clusters than Li
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
                               [&](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
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
                               {"MergerData", "TPCData"});

    auto dfZandQ = dfFilter.Define("zDrift",
                                   [&](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                                   {
                                       auto& rp {m.fRP};

                                       auto idx = m.fLightIdx;
                                       if(idx < 0)
                                           return -100.0f;

                                       if(static_cast<std::size_t>(idx) >= tpc.fClusters.size())
                                           return -100.0f;

                                       auto& cluster = tpc.fClusters.at(idx);
                                       auto& voxels = cluster.GetRefToVoxels();

                                       if(tpc.fRPs.empty())
                                           return -100.0f;
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

                                       float zDrift = 110.0 + deltaZ;

                                       return zDrift;
                                   },
                                   {"MergerData", "TPCData"});

    auto dfZandQAux = dfFilterAux.Define("zDrift",
                                         [&](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                                         {
                                             auto& rp {m.fRP};

                                             auto idx = m.fLightIdx;
                                             if(idx < 0)
                                                 return -100.0f;

                                             if(static_cast<std::size_t>(idx) >= tpc.fClusters.size())
                                                 return -100.0f;

                                             auto& cluster = tpc.fClusters.at(idx);
                                             auto& voxels = cluster.GetRefToVoxels();

                                             if(tpc.fRPs.empty())
                                                 return -100.0f;
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

                                             float zDrift = 110.0 + deltaZ;

                                             return zDrift;
                                         },
                                         {"MergerData", "TPCData"});

    struct ZDriftCut
    {
        double zMin; // meanStart + nSigma * sigmaStart
        double zMax; // meanEnd   - nSigma * sigmaEnd
    };
    const double nSigmaZDrift {3.0};
    std::map<int, ZDriftCut> zDriftCuts;
    {
        std::ifstream finZ("./Outputs/zDrift_cut_perRun.dat");
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

    auto dfFilterZandQ = dfZandQ.Filter(
        [&zDriftCuts](float zDrift, ActRoot::MergerData& m)
        {
            auto it = zDriftCuts.find(m.fRun);
            if(it == zDriftCuts.end())
                return false; // run not in file → reject
            return zDrift > it->second.zMin && zDrift < it->second.zMax && m.fLight.fQave > 600;
        },
        {"zDrift", "MergerData"});

    auto dfFilterZandQAux = dfZandQAux.Filter(
        [&zDriftCuts](float zDrift, ActRoot::MergerData& m)
        {
            auto it = zDriftCuts.find(m.fRun);
            if(it == zDriftCuts.end())
                return false; // run not in file → reject
            return zDrift > it->second.zMin && zDrift < it->second.zMax && m.fLight.fQave > 600;
        },
        {"zDrift", "MergerData"});

    ActRoot::CutsManager<std::string> cuts;
    cuts.ReadCut("l1", TString::Format("./Cuts/%s_TLvsQ_%s.root", light.c_str(), beam.c_str()).Data());
    cuts.ReadCut("deuterium", "./Cuts/elasticInproton_Theta_Q.root");
    auto dfElasticInProtonBanana = dfZandQ.Filter(
        [&](ActRoot::MergerData& m)
        {
            return cuts.IsInside("l1", m.fLight.fRawTL, m.fLight.fQtotal) &&
                   cuts.IsInside("deuterium", m.fThetaLight, m.fLight.fQtotal);
        },
        {"MergerData"});

    auto dfElasticInProtonBananaAux = dfZandQAux.Filter(
        [&](ActRoot::MergerData& m)
        {
            return cuts.IsInside("l1", m.fLight.fRawTL, m.fLight.fQtotal) &&
                   cuts.IsInside("deuterium", m.fThetaLight, m.fLight.fQtotal);
        },
        {"MergerData"});

    std::cout << "Events in proton banana rebined: " << *dfElasticInProtonBanana.Count() << "\n";
    std::cout << "Events in proton banana: " << *dfElasticInProtonBananaAux.Count() << "\n";

    // Check number of entries with differences between dfElasticInProtonBanana and dfElasticInProtonBananaAux comparing TL and Qtotal
    





    
    // Plot TL and Qtotal for each df to check if something changed
    auto hTL = dfElasticInProtonBanana.Histo1D(
        {"hTL", "Track Length;Track Length [a.u.]", 240, 0, 120},
        "MergerData.fLight.fRawTL");
    auto hTL_Aux = dfElasticInProtonBananaAux.Histo1D(
        {"hTL_Aux", "Track Length;Track Length [a.u.]", 240, 0, 120},
        "MergerData.fLight.fRawTL");
    auto hQtot = dfElasticInProtonBanana.Histo1D(
        {"hQtot", "Qtotal;Qtotal [a.u.]", 200, 0, 3e5},
        "MergerData.fLight.fQtotal");
    auto hQtot_Aux = dfElasticInProtonBananaAux.Histo1D(
        {"hQtot_Aux", "Qtotal;Qtotal [a.u.]", 200, 0, 3e5},
        "MergerData.fLight.fQtotal");

    auto* c1 = new TCanvas("c1", "c1", 1200, 600);
    c1->Divide(2, 2);
    c1->cd(1);
    hTL->DrawClone("COLZ");
    c1->cd(2);
    hTL_Aux->DrawClone("COLZ");
    c1->cd(3);
    hQtot->DrawClone("COLZ");
    c1->cd(4);
    hQtot_Aux->DrawClone("COLZ");
}