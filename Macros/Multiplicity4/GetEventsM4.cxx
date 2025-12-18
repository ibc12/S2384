#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActSRIM.h"
#include "ActSilData.h"
#include "ActSilMatrix.h"
#include "ActSilSpecs.h"
#include "ActTPCData.h"
#include "ActParticle.h"
#include "ActKinematics.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"

#include <fstream>

float scaleXY = 2.0f;    // mm per bin
float scaleZ = 2.84032f; // mm per bin

void ScalePoint(ROOT::Math::XYZPointF& point, bool addOffset = false)
{
    if(addOffset) // when converting a bin point to physical units which wasnt already corrected
        point += ROOT::Math::XYZVector {0.5, 0.5, 0.5};
    point.SetX(point.X() * scaleXY);
    point.SetY(point.Y() * scaleXY);
    point.SetZ(point.Z() * scaleZ);
}


void GetEventsM4()
{
    // As I did it before I depend on the Merger to get the data, with that I cannot recover the f0 information, because
    // the merger erase events with multiplicity two in the slicons

    std::string dataconf = "../../configs/data_7Li.conf";

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

    df.Describe().Print();

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
                            {"TPCData", "GETTree.TPCData"});

    // Get silicon matrix and silspecs to get particle to go to lateral silicons
    std::string filenameSMleft {"../SilVetos/Outputs/Dists/sms_l0.root"};
    auto fileSMleft {new TFile {filenameSMleft.c_str()}};
    ActPhysics::SilMatrix* smleft =
        fileSMleft->Get<ActPhysics::SilMatrix>("sm5"); // matrix for good distance of left wall
    double silCentreLeft = smleft->GetMeanZ({4, 5});
    double beamOffsetLeft {3.36}; // mm offset of beam with respect to sils 4 and 5 off left wall (need to lower beam
                                  // that amount respect of silicons)
    const double zVertexMeanLeft {silCentreLeft - beamOffsetLeft}; // beam on left
    std::string filenameSMright {"../SilVetos/Outputs/Dists/sms_r0.root"};
    auto fileSMright {new TFile {filenameSMright.c_str()}};
    ActPhysics::SilMatrix* smright =
        fileSMright->Get<ActPhysics::SilMatrix>("sm5"); // matrix for good distance of right wall
    double silCentreRight = smright->GetMeanZ({4, 5});
    double beamOffsetRight {-2.14};
    const double zVertexMeanRight {silCentreRight - beamOffsetRight}; // beam on right

    auto sils {std::make_shared<ActPhysics::SilSpecs>()};
    sils->ReadFile("../../configs/silspecs.conf");


    // Now get the index of the particle that hit the silicon to do PID
    auto dfWithIndex = dfFilter.Define("Layer", [](ActRoot::ModularData& mod) { 
        int gatConf = mod.Get("GATCONF");
        std::string layer {""};
        if(gatConf == 1)
            layer = "l0";
        else if(gatConf == 2)
            layer = "r0";
        else if(gatConf == 4)
            layer = "f0";
        return layer;
    }, {"ModularData"}).Define("pseudoSP", [&](ActRoot::TPCData& tpc, ActRoot::SilData& sil, ActRoot::ModularData& mod, std::string layer){
        auto pseudoSP = ROOT::Math::XYZPointF {};
        if(layer == "")
            return pseudoSP;

        sil.ApplyFinerThresholds(sils);
        // Now get aprox point of silicon
        auto pointLayer {sils->GetLayer(layer).GetPoint()};
        // Move XZ coordinates to correct sil unit
        auto silIdx {sil.fSiN.at(layer)[0]};
        auto silHitPosition {sils->GetLayer(layer).GetPlacements().at(silIdx)};
     
        if(layer == "l0")
            pseudoSP = {static_cast<float>(silHitPosition.first), static_cast<float>(pointLayer.Y()),
                        static_cast<float>(smleft->GetMeanZ({silIdx}) - beamOffsetLeft)}; // this is all in mm
        else if(layer == "r0")  // r0
        {
            pseudoSP = {static_cast<float>(silHitPosition.first), static_cast<float>(pointLayer.Y()),
                        static_cast<float>(smright->GetMeanZ({silIdx}) - beamOffsetRight)}; // this is all in mm
        }
        return pseudoSP;
    }, {"TPCData", "SilData", "ModularData", "Layer"})
    .Define("ParticleIndex", [&](ActRoot::TPCData& tpc, ActRoot::SilData& sil, ActRoot::ModularData& mod, std::string layer, ROOT::Math::XYZPointF pseudoSP)
                                       { // Get index wall that triggered
        if(layer == "")
            return -1;

        sil.ApplyFinerThresholds(sils);
        
            // Now move each cluster line to pseudoSP x and get distance, lower distance means it hit
            float minDistance = 1e8;
            int bestIdx = -1;
            for(int i = 0; i < tpc.fClusters.size(); ++i)
            {
                auto& cluster = tpc.fClusters[i];

                auto line = cluster.GetLine();
                line.Scale(scaleXY, scaleZ);

                auto posibleSP = line.MoveToY(pseudoSP.Y());
                auto distance = (posibleSP - pseudoSP).R();
                // std::cout << "Cluster " << i << " distance to silicon hit point: " << distance << " mm" << std::endl;

                if(distance < minDistance)
                {
                    minDistance = distance;
                    bestIdx = i;
                }
            }
            return bestIdx;
        }, {"TPCData", "SilData", "ModularData", "Layer"}).Define("DistanceInGas", [&](ActRoot::TPCData& tpc, int idx, ROOT::Math::XYZPointF pseudoSP)
                                       { // Get distance in gas to that wall
        if(pseudoSP == ROOT::Math::XYZPointF{}) // Check if pseudoSP is default constructed
            return -1.f;
        auto& cluster = tpc.fClusters[idx];
        auto line = cluster.GetLine();
        line.Scale(scaleXY, scaleZ);
        auto rp = tpc.fRPs.front();
        auto sp = line.MoveToY(pseudoSP.Y());
        auto distance = (sp - pseudoSP).R();
        return distance;
        }, {"TPCData", "ParticleIndex", "pseudoSP"});

    auto hIdx = dfWithIndex.Histo1D({"hIdx", "Index of particle hitting lateral silicon;Index;Counts", 7, -2, 5},
                                    "ParticleIndex");
    TCanvas* c0 = new TCanvas("c0", "Index of particle hitting lateral silicon", 800, 600);
    hIdx->DrawClone();
    // Now define Qave for that index to do PID plot
    auto dfWithQave = dfWithIndex
                          .Define("QaveSelectedParticle",
                                  [](ActRoot::TPCData& tpc, int idx)
                                  {
                                      auto cluster = tpc.fClusters[idx];
                                      cluster.SortAlongDir(cluster.GetLine().GetDirection());
                                      auto voxels = cluster.GetRefToVoxels();
                                      float Q = 0.f;
                                      for(const auto& v : voxels)
                                          Q += v.GetCharge();
                                      ROOT::Math::XYZPointF firstPos = voxels.front().GetPosition();
                                      ROOT::Math::XYZPointF lastPos = voxels.back().GetPosition();
                                      ScalePoint(firstPos);
                                      ScalePoint(lastPos);
                                      float length = (lastPos - firstPos).R();
                                      if(length <= 0.0)
                                          return -1.f;
                                      return Q / length;
                                  },
                                  {"TPCData", "ParticleIndex"})
                          .Define("SilESelectedParticle",
                                  [&](ActRoot::SilData& sil, ActRoot::ModularData& mod, std::string layer)
                                  {
                                      std::cout << "------------------------" << std::endl;
                                      std::cout << "Size of hits in layer  before threshold" << layer << ": "
                                                << sil.fSiE.at(layer).size() << std::endl;
                                      sil.ApplyFinerThresholds(sils);
                                      std::cout << "After threshold: " << sil.fSiE.at(layer).size() << std::endl;
                                      return sil.fSiE.at(layer)[0];
                                  },
                                  {"SilData", "ModularData", "Layer"});

    // Plot PID to check everything is ok, first divide in sil layers
    auto dfL0 = dfWithQave.Filter([](ActRoot::ModularData& mod) { return mod.Get("GATCONF") == 1; }, {"ModularData"});
    auto dfR0 = dfWithQave.Filter([](ActRoot::ModularData& mod) { return mod.Get("GATCONF") == 2; }, {"ModularData"});

    auto hPIDl0 =
        dfL0.Histo2D({"hPIDl0", "PID plot Qave vs SilE;Silicon Energy (MeV);Q_{ave} (a.u.)", 100, 0, 30, 100, 0, 500},
                     "SilESelectedParticle", "QaveSelectedParticle");
    auto hPIDr0 =
        dfR0.Histo2D({"hPIDr0", "PID plot Qave vs SilE;Silicon Energy (MeV);Q_{ave} (a.u.)", 100, 0, 30, 100, 0, 500},
                     "SilESelectedParticle", "QaveSelectedParticle");
    // Get also the cuts for p and d, and see if events are inside them
    auto cuts = ActRoot::CutsManager<std::string> {};
    cuts.ReadCut("l0_p", TString::Format("../../PostAnalysis/Cuts/pid_p_l0_7Li.root").Data());
    cuts.ReadCut("r0_p", TString::Format("../../PostAnalysis/Cuts/pid_p_r0_7Li.root").Data());
    cuts.ReadCut("l0_d", TString::Format("../../PostAnalysis/Cuts/pid_d_l0_7Li.root").Data());
    cuts.ReadCut("r0_d", TString::Format("../../PostAnalysis/Cuts/pid_d_r0_7Li.root").Data());

    TCanvas* c1 = new TCanvas("c1", "PID plot", 1600, 1200);
    c1->Divide(2, 1);
    c1->cd(1);
    hPIDl0->DrawClone("colz");
    cuts.DrawCut("l0_p");
    cuts.DrawCut("l0_d");
    c1->cd(2);
    hPIDr0->DrawClone("colz");
    cuts.DrawCut("r0_p");
    cuts.DrawCut("r0_d");

    std::cout << "Number of events with M=4 and later silicon hit: " << dfFilter.Count().GetValue() << std::endl;

    // Reconstruct ex with energy of protons or deuterons at vertex
    // For that we need the srim files n the gas for the light particles
    auto* srim {new ActPhysics::SRIM};
    srim->ReadTable("p", TString::Format("../Calibrations/SRIM/1H_900mb_CF4_95-5.txt").Data());
    srim->ReadTable("d", TString::Format("../Calibrations/SRIM/2H_900mb_CF4_95-5.txt").Data());
    srim->ReadTable("7Li", TString::Format("../Calibrations/SRIM/7Li_900mb_CF4_95-5.txt").Data());
    srim->ReadTable("mylar", TString::Format("../Calibrations/SRIM/7Li_Mylar.txt").Data());

    ActPhysics::Particle pBeam {"7Li"};
    ActPhysics::Particle pDeuteron {"d"};
    ActPhysics::Particle pProton {"p"};

     // Initial energy
    double initialEnergy {7.558}; // meassured by operators; resolution of 0,19%
    initialEnergy = srim->Slow("mylar", initialEnergy * pBeam.GetAMU(), 0.0168);
    initialEnergy = srim->Slow("7Li", initialEnergy, 60); // 60 mm of gas before the pad plane
    initialEnergy = initialEnergy / pBeam.GetAMU(); // back to amu units

    auto dfVertex = dfWithQave.Define("EVertex",
                                      [&](float distGas, float eSil)
                                      {
                                          double ret = srim->EvalInitialEnergy("p", eSil, distGas);
                                          // else // L1 trigger
                                          //     ret = srim->EvalEnergy("p", d.fLight.fTL);
                                          return ret;
                                      },
                                      {"DistanceInGas", "SilESelectedParticle"});

    ActPhysics::Kinematics kinP {pBeam, pDeuteron, pProton, initialEnergy * pBeam.GetAMU()};
    ActPhysics::Kinematics kinD {pBeam, pProton, pDeuteron, initialEnergy * pBeam.GetAMU()};
    std::vector<ActPhysics::Kinematics> vkins {dfVertex.GetNSlots()};
    for(auto& k : vkins)
        k = kin;
    dfVertex =
        dfVertex.DefineSlot("Ex",
                       [&](unsigned int slot, const ActRoot::MergerData& d, double EVertex, double EBeam)
                       {
                           vkins[slot].SetBeamEnergy(EBeam);
                           return vkins[slot].ReconstructExcitationEnergy(EVertex, (d.fThetaLight) * TMath::DegToRad());
                       },
                       {"MergerData", "EVertex", "EBeam"});
    dfVertex =
        dfVertex.DefineSlot("ThetaCM",
                       [&](unsigned int slot, const ActRoot::MergerData& d, double EVertex, double EBeam)
                       {
                           vkins[slot].SetBeamEnergy(EBeam);
                           return vkins[slot].ReconstructTheta3CMFromLab(EVertex, (d.fThetaLight) * TMath::DegToRad()) *
                                  TMath::RadToDeg();
                       },
                       {"MergerData", "EVertex", "EBeam"});

    // Forech to the events with output
    // std::ofstream out("./Outputs/PIDoutCuts.dat");
    // dfWithQave.Foreach(
    //     [&](ActRoot::MergerData& m, float silE)
    //     {
    //         if(silE < 0.6)
    //             m.Stream(out);
    //     },
    //     {"MergerData", "SilESelectedParticle"});
    // out.close();
}