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

void CheckHitsLatNoLatPID()
{
    // Get Output from pipe0
    ROOT::RDataFrame df("Final_Tree", TString::Format("../Outputs/SelectorM4_7Li.root").Data());

    auto sils {std::make_shared<ActPhysics::SilSpecs>()};
    sils->ReadFile("../../../configs/silspecs.conf");

    // Get silicon matrix and silspecs to get particle to go to lateral silicons
    std::string filenameSMleft {"../../SilVetos/Outputs/Dists/sms_l0.root"};
    auto fileSMleft {new TFile {filenameSMleft.c_str()}};
    ActPhysics::SilMatrix* smleft =
        fileSMleft->Get<ActPhysics::SilMatrix>("sm5"); // matrix for good distance of left wall
    double silCentreLeft = smleft->GetMeanZ({4, 5});
    double beamOffsetLeft {3.36}; // mm offset of beam with respect to sils 4 and 5 off left wall (need to lower beam
                                  // that amount respect of silicons)
    const double zVertexMeanLeft {silCentreLeft - beamOffsetLeft}; // beam on left
    std::string filenameSMright {"../../SilVetos/Outputs/Dists/sms_r0.root"};
    auto fileSMright {new TFile {filenameSMright.c_str()}};
    ActPhysics::SilMatrix* smright =
        fileSMright->Get<ActPhysics::SilMatrix>("sm5"); // matrix for good distance of right wall
    double silCentreRight = smright->GetMeanZ({4, 5});
    double beamOffsetRight {-2.14};
    const double zVertexMeanRight {silCentreRight - beamOffsetRight}; // beam on right
    std::string filenameSMfront {"../../SilVetos/Outputs/Dists/sms_f0.root"};
    auto fileSMfront {new TFile {filenameSMfront.c_str()}};
    ActPhysics::SilMatrix* smfront =
        fileSMfront->Get<ActPhysics::SilMatrix>("sm2"); // matrix for good distance of front wall
    double silCentreFront = smfront->GetMeanZ({7, 6, 4});
    double beamOffsetFront {-2.5};
    const double zVertexMeanFront {silCentreFront - beamOffsetFront}; // beam on front

    std::cout << "Silicon left centre at Z = " << silCentreLeft << " mm" << std::endl;
    std::cout << "Silicon right centre at Z = " << silCentreRight << " mm" << std::endl;
    std::cout << "Silicon front centre at Z = " << silCentreFront << " mm" << std::endl;

    // First Get event 18512 r66
    df.Filter([](ActRoot::MergerData m) { return m.fRun == 66 && m.fEntry == 20187; }, {"MergerData"})
        .Foreach(
            [&](ActRoot::MergerData m, ActRoot::TPCData tpc, int LightIdx, ActRoot::SilData silIn)
            {
                std::cout << "The Light Cluster is: " << LightIdx << std::endl;
                // Print the angle and the endpoint of all clusters
                for(int i = 0; i < tpc.fClusters.size(); ++i)
                {
                    auto& cluster = tpc.fClusters[i];
                    auto line = cluster.GetLine();
                    auto dir = line.GetDirection();
                    cluster.SortAlongDir(dir);
                    auto voxels = cluster.GetRefToVoxels();
                    auto lastPos = voxels.back().GetPosition();
                    auto theta = Utils::GetTheta3D({1, 0, 0}, dir);
                    std::cout << "-------------------------" << std::endl;
                    std::cout << "Cluster " << i << ": theta = " << theta << " degrees" << std::endl;
                    std::cout << "Cluster " << i << ": endpoint = (" << lastPos.X() << ", " << lastPos.Y() << ", "
                              << lastPos.Z() << ")" << std::endl;
                    std::cout << " --------------------------" << std::endl;
                }

                // Get distance to hit and see why fails
                ActRoot::SilData sil {silIn};
                sil.ApplyFinerThresholds(sils);

                auto& cluster = tpc.fClusters[LightIdx];
                auto line = cluster.GetLine();
                line.Scale(Utils::scaleXY, Utils::scaleZ);

                double minDistance {1e8};
                const double kSilMatchTol {
                    60.}; // mm -- distance between front and sides is 71 mm. So has to be lower than that

                std::cout << "Distance to hit in lateral silicons: " << std::endl;

                for(const std::string layer : {"l0", "r0", "f0"})
                {
                    std::cout << " --------------------------" << std::endl;
                    std::cout << "Checking layer: " << layer << std::endl;
                    auto it = sil.fSiN.find(layer);
                    if(it == sil.fSiN.end() || it->second.empty())
                        continue; // no hit en esta layer

                    auto pointLayer {sils->GetLayer(layer).GetPoint()};
                    auto silIdx {it->second[0]};

                    auto silHitPosition {sils->GetLayer(layer).GetPlacements().at(silIdx)};
                    std::cout << "Silicon idx: " << silIdx << " hit position: (" << silHitPosition.first << ", "
                              << silHitPosition.second << std::endl;

                    ROOT::Math::XYZPointF pseudoSP {};
                    if(layer == "l0")
                        pseudoSP = {static_cast<float>(silHitPosition.first), static_cast<float>(pointLayer.Y()),
                                    static_cast<float>(smleft->GetMeanZ({silIdx}) - beamOffsetLeft)};
                    else if(layer == "r0")
                        pseudoSP = {static_cast<float>(silHitPosition.first), static_cast<float>(pointLayer.Y()),
                                    static_cast<float>(smright->GetMeanZ({silIdx}) - beamOffsetRight)};
                    else if(layer == "f0")
                        pseudoSP = {static_cast<float>(silHitPosition.first), static_cast<float>(pointLayer.Y()),
                                    static_cast<float>(smfront->GetMeanZ({silIdx}) - beamOffsetFront)};
                    // if(layer == "l0")
                    //     pseudoSP = {static_cast<float>(pointLayer.X()), static_cast<float>(pointLayer.Y()), 0};
                    // else if(layer == "r0")
                    //     pseudoSP = {static_cast<float>(pointLayer.X()), static_cast<float>(pointLayer.Y()), 0};
                    // else if(layer == "f0")
                    //     pseudoSP = {static_cast<float>(pointLayer.X()), static_cast<float>(pointLayer.Y()), 0};
                    std::cout << "PseudoSP for layer " << layer << ": (" << pseudoSP.X() << ", " << pseudoSP.Y() << ", "
                              << pseudoSP.Z() << ")" << std::endl;

                    ROOT::Math::XYZPoint posibleSP = {0, 0, 0};
                    float distance {};
                    // if(layer == "l0" || layer == "r0")
                    //     distance = (posibleSP - pseudoSP).R();
                    if(layer == "l0" || layer == "r0")
                    {
                        posibleSP = line.MoveToY(pointLayer.Y());
                        // posibleSP.SetZ(0); // We have a bad matrix for the front silicon, so just check in XY plane
                        distance = (posibleSP - pseudoSP).R();
                    }
                    else if(layer == "f0")
                    {
                        posibleSP = line.MoveToX(pointLayer.X());
                        // posibleSP.SetZ(0); // We have a bad matrix for the front silicon, so just check in XY plane
                        distance = (posibleSP - pseudoSP).R();
                    }
                    std::cout << "PosibleSP for layer " << layer << ": (" << posibleSP.X() << ", " << posibleSP.Y()
                              << ", " << posibleSP.Z() << ")" << std::endl;

                    std::cout << "Distance to layer " << layer << ": " << distance << " mm" << std::endl;
                    std::cout << " --------------------------" << std::endl;
                }
            },
            {"MergerData", "TPCData", "LightIdx", "SilData"});
}