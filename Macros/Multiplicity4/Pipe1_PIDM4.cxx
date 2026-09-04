#ifndef Pipe1_PIDM4_cxx
#define Pipe1_PIDM4_cxx

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
#include "TMath.h"
#include "TString.h"

#include <map>
#include <string>

#include "./Utils.h"

// Result of geometrically matching the light-particle track (chosen by angle,
// see LightIdx below) to a real hit in either lateral silicon wall, regardless
// of which wall actually fired the trigger.
struct LightMatch
{
    std::string layer {""};
    ROOT::Math::XYZPointF pseudoSP {};
    float distance {1e8};
};

void Pipe1_PIDM4(const std::string& beam, const std::string& target, const std::string& light)
{
    // Get file from pipe0
    TString infile = TString::Format("./Outputs/SelectorM4_%s.root", beam.c_str());
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df("Final_Tree", infile.Data());

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


    // Silicon side length is 50 mm; use half of that as the matching tolerance
    // between a track extrapolated to a silicon plane and the actual hit point.
    const double kSilMatchTolerance {100.}; // mm

    // Step 1: identify the beam-like cluster (independent of trigger/layer).
    auto dfWithIndex =
        df.Define("BeamIdx",
                  [](ActRoot::TPCData& tpc)
                  { // Get index of beam-like particle
                      int beamIdx = -1;
                      for(int i = 0; i < tpc.fClusters.size(); ++i)
                      {
                          auto& cluster = tpc.fClusters[i];
                          if(cluster.GetIsBeamLike())
                          {
                              beamIdx = i;
                              break;
                          }
                      }
                      return beamIdx;
                  },
                  {"TPCData"})
            // Step 2: the light particle (p, d, ...) is, almost always, the
            // non-beam cluster that makes the largest angle with the beam
            // direction. This no longer depends on which silicon fired the
            // trigger, so it also works for frontal-trigger (GATCONF==4)
            // events where the light particle still reaches a lateral wall.
            .Define("LightIdx",
                    [](ActRoot::TPCData& tpc, int beamIdx)
                    {
                        if(beamIdx < 0)
                            return -1;
                        auto beamDir = tpc.fClusters[beamIdx].GetLine().GetDirection().Unit();
                        double maxAngle {-1.};
                        int bestIdx {-1};
                        for(int i = 0; i < tpc.fClusters.size(); ++i)
                        {
                            if(i == beamIdx)
                                continue;
                            auto dir = tpc.fClusters[i].GetLine().GetDirection().Unit();
                            double angle = TMath::ACos(beamDir.Dot(dir)) * TMath::RadToDeg();
                            if(angle > maxAngle)
                            {
                                maxAngle = angle;
                                bestIdx = i;
                            }
                        }
                        return bestIdx;
                    },
                    {"TPCData", "BeamIdx"})
            // Step 3: check, purely geometrically, whether the extrapolated
            // LightIdx track matches a real hit in EITHER lateral silicon wall
            // (l0 or r0), independently of GATCONF. Keep whichever wall gives
            // the smallest extrapolation-to-hit distance, provided it is
            // within one silicon half-width (25 mm) of the real hit.
            .Define("LightMatch",
                    [&](ActRoot::TPCData& tpc, ActRoot::SilData& sil, int lightIdx) -> LightMatch
                    {
                        LightMatch match {};
                        if(lightIdx < 0 || lightIdx >= (int)tpc.fClusters.size())
                            return match;

                        sil.ApplyFinerThresholds(sils);

                        auto& cluster = tpc.fClusters[lightIdx];
                        auto line = cluster.GetLine();
                        line.Scale(Utils::scaleXY, Utils::scaleZ);

                        for(const auto& layer : {std::string {"l0"}, std::string {"r0"}})
                        {
                            if(!sil.fSiN.count(layer) || sil.fSiN.at(layer).empty())
                                continue; // no hit registered in this wall for this event

                            auto pointLayer {sils->GetLayer(layer).GetPoint()};
                            auto silIdx {sil.fSiN.at(layer)[0]};
                            auto silHitPosition {sils->GetLayer(layer).GetPlacements().at(silIdx)};

                            ROOT::Math::XYZPointF pseudoSP {};
                            if(layer == "l0")
                                pseudoSP = {static_cast<float>(silHitPosition.first),
                                            static_cast<float>(pointLayer.Y()),
                                            static_cast<float>(smleft->GetMeanZ({silIdx}) - beamOffsetLeft)};
                            else // r0
                                pseudoSP = {static_cast<float>(silHitPosition.first),
                                            static_cast<float>(pointLayer.Y()),
                                            static_cast<float>(smright->GetMeanZ({silIdx}) - beamOffsetRight)};

                            auto posibleSP = line.MoveToY(pseudoSP.Y());
                            auto distance = (posibleSP - pseudoSP).R();

                            if(distance < match.distance)
                            {
                                match.distance = distance;
                                match.layer = layer;
                                match.pseudoSP = pseudoSP;
                            }
                        }

                        if(match.distance > kSilMatchTolerance)
                        {
                            match.layer = "";
                            match.pseudoSP = ROOT::Math::XYZPointF {};
                        }

                        return match;
                    },
                    {"TPCData", "SilData", "LightIdx"})
            .Define("Layer", "LightMatch.layer")
            .Define("pseudoSP", "LightMatch.pseudoSP")
            // Drop events where the light particle's track doesn't geometrically
            // match a hit in either lateral wall (no reliable silicon PID/energy
            // is possible for these).
            .Filter([](const std::string& layer) { return layer != ""; }, {"Layer"})
            .Define("TrackLength",
                    [&](ActRoot::TPCData& tpc, int idx, ROOT::Math::XYZPointF pseudoSP)
                    { // Get distance in gas to that wall
                        auto& cluster = tpc.fClusters[idx];
                        auto line = cluster.GetLine();
                        line.Scale(Utils::scaleXY, Utils::scaleZ);
                        auto rp = tpc.fRPs.front();
                        auto sp = line.MoveToY(pseudoSP.Y());
                        auto distance = (sp - rp).R();
                        return distance;
                    },
                    {"TPCData", "LightIdx", "pseudoSP"});

    auto hIdx = dfWithIndex.Histo1D(
        {"hIdx", "Idx of max-angle (light) particle, matched to a lateral silicon;Index;Counts", 7, -2, 5},
        "LightIdx");

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
                                      Utils::ScalePoint(firstPos);
                                      Utils::ScalePoint(lastPos);
                                      float length = (lastPos - firstPos).R();
                                      if(length <= 0.0)
                                          return -1.f;
                                      return Q / length;
                                  },
                                  {"TPCData", "LightIdx"})
                          .Define("SilESelectedParticle",
                                  [&](ActRoot::SilData& sil, ActRoot::ModularData& mod, std::string layer)
                                  {
                                      // std::cout << "------------------------" << std::endl;
                                      // std::cout << "Size of hits in layer  before threshold" << layer << ": "
                                      //           << sil.fSiE.at(layer).size() << std::endl;
                                      sil.ApplyFinerThresholds(sils);
                                      // std::cout << "After threshold: " << sil.fSiE.at(layer).size() << std::endl;
                                      return sil.fSiE.at(layer)[0];
                                  },
                                  {"SilData", "ModularData", "Layer"})
                          .Define("ThetaLab",
                                  [&](ActRoot::TPCData& tpc, int Lightidx, int BeamIdx)
                                  {
                                      auto cluster = tpc.fClusters[Lightidx];
                                      auto beam = tpc.fClusters[BeamIdx];
                                      auto theta =
                                          Utils::GetTheta3D(beam.GetLine().GetDirection(), cluster.GetLine().GetDirection());
                                      return theta; // in degrees
                                  },
                                  {"TPCData", "LightIdx", "BeamIdx"});

    // Plot PID to check everything is ok, first divide by the geometrically
    // matched lateral wall (Layer), not by which wall triggered.
    auto dfL0 = dfWithQave.Filter([](const std::string& layer) { return layer == "l0"; }, {"Layer"});
    auto dfR0 = dfWithQave.Filter([](const std::string& layer) { return layer == "r0"; }, {"Layer"});

    auto hPIDl0 =
        dfL0.Histo2D({"hPIDl0", "PID plot Qave vs SilE;Silicon Energy (MeV);Q_{ave} (a.u.)", 100, 0, 30, 100, 0, 500},
                     "SilESelectedParticle", "QaveSelectedParticle");
    auto hPIDr0 =
        dfR0.Histo2D({"hPIDr0", "PID plot Qave vs SilE;Silicon Energy (MeV);Q_{ave} (a.u.)", 100, 0, 30, 100, 0, 500},
                     "SilESelectedParticle", "QaveSelectedParticle");

    auto cuts = ActRoot::CutsManager<std::string> {};
    cuts.ReadCut("l0",
                 TString::Format("../../PostAnalysis/Cuts/pid_%s_l0_%s.root", light.c_str(), beam.c_str()).Data());
    cuts.ReadCut("r0",
                 TString::Format("../../PostAnalysis/Cuts/pid_%s_r0_%s.root", light.c_str(), beam.c_str()).Data());

    std::cout << "Number of events with M=4 and later silicon hit: " << dfWithQave.Count().GetValue() << std::endl;

    // Finally filter the events with the cuts
    auto dfPIDfiltered = dfWithQave.Filter([&cuts](std::string layer, float silE, float Qave)
                                   {
                                       if(layer == "l0")
                                           return cuts.IsInside("l0", silE, Qave);
                                       else
                                           return cuts.IsInside("r0", silE, Qave);
                                   },
                                   {"Layer", "SilESelectedParticle", "QaveSelectedParticle"});

    TCanvas* c0 = new TCanvas("Pipe1_0", "Index of particle hitting lateral silicon", 800, 600);
    hIdx->DrawClone();

    TCanvas* c1 = new TCanvas("Pipe1_1", "PID plot", 1600, 1200);
    c1->Divide(2, 1);
    c1->cd(1);
    hPIDl0->DrawClone("colz");
    cuts.DrawCut("l0");
    c1->cd(2);
    hPIDr0->DrawClone("colz");
    cuts.DrawCut("r0");

    // Save dataframe in a .root file
    TString outfile = TString::Format("./Outputs/PIDM4_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str());
    dfPIDfiltered.Snapshot("Final_Tree", outfile);
    std::cout << "Saving Final_Tree in " << outfile << " from Pipe1_PIDM4" << '\n';
}
#endif