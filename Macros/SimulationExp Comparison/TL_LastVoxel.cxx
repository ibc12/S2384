#include "ActMergerData.h"
#include "ActTPCData.h"
#include "ActInputParser.h"

#include "ROOT/RDataFrame.hxx"

#include "TString.h"
#include "TCanvas.h"

#include <string>
#include <vector>


void TL_LastVoxel()
{
    // Get data for 7Li (d,d) elastic scattering
    std::string beam {"7Li"};
    std::string target {"2H"};
    std::string targetExp {"d"};
    std::string light {"2H"};
    std::string lightExp {"d"};
    std::string heavy {"7Li"};
    double Ex {0.}; // MeV
    // Get simu data
    TString fileName {TString::Format("../../Simulation/Outputs/%s/%s_%s_TRIUMF_Eex_%.3f_nPS_%d_pPS_%d_L1.root",
                                      beam.c_str(), target.c_str(), light.c_str(), Ex, 0, 0)};
    // Get exp data
    TString fileNameExp {TString::Format("../../PostAnalysis/Outputs/tree_ex_F_%s_%s_%s.root", beam.c_str(),
                                         targetExp.c_str(), lightExp.c_str())};

    auto df = ROOT::RDataFrame("SimulationTTree", fileName);
    auto dfExp = ROOT::RDataFrame("Final_Tree", fileNameExp);

    // Filter exp data to get only L1 events
    auto dfExpL1 = dfExp.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"});

    // Get drift parameter from config file
    ActRoot::InputParser parser {};
    parser.ReadFile("../../configs/detector.conf");
    auto driftBlock = parser.GetBlock("Merger");
    auto driftFactor = driftBlock->GetDouble("DriftFactor"); // in mm^2/us

    // Define for exp data the position x and y of the last voxel
    auto dfExpL1Positions = dfExpL1
                                .Define("LastPosX",
                                        [&driftFactor](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                                        {
                                            auto lightIdx = m.fLightIdx;
                                            auto lightCl = tpc.fClusters[lightIdx];
                                            auto line = lightCl.GetLine();
                                            auto voxels = lightCl.GetVoxels();
                                            auto vertex = tpc.fRPs.front();
                                            auto gravityPoint = line.GetPoint();
                                            auto directionToSort = gravityPoint - vertex;
                                            // Sort voxels according to their projection on the line
                                            lightCl.SortAlongDir(directionToSort);
                                            lightCl.ReFit();
                                            auto lastVoxel = lightCl.GetVoxels().back();
                                            auto lastVoxelPos = lastVoxel.GetPosition();
                                            return lastVoxelPos.X() * 2; // Return in mm
                                        },
                                        {"MergerData", "TPCData"})
                                .Define("LastPosY",
                                        [&driftFactor](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                                        {
                                            auto lightIdx = m.fLightIdx;
                                            auto lightCl = tpc.fClusters[lightIdx];
                                            auto line = lightCl.GetLine();
                                            auto voxels = lightCl.GetVoxels();
                                            auto vertex = tpc.fRPs.front();
                                            auto gravityPoint = line.GetPoint();
                                            auto directionToSort = gravityPoint - vertex;
                                            // Sort voxels according to their projection on the line
                                            lightCl.SortAlongDir(directionToSort);
                                            lightCl.ReFit();
                                            auto lastVoxel = lightCl.GetVoxels().back();
                                            auto lastVoxelPos = lastVoxel.GetPosition();
                                            return lastVoxelPos.Y() * 2; // Return in mm
                                        },
                                        {"MergerData", "TPCData"});


    // Get the TL distribution for simu
    auto hTLSimu = df.Histo1D({"hSimu", "TL distribution for simu;TL (mm);Counts", 200, 0, 200}, "TL");
    auto hLastVoxelSimu = df.Histo2D(
        {"hLastVoxelSimu", "Last voxel position for simu;LastPosX (mm);LastPosY (mm)", 256 / 2, 0, 256, 256 / 2, 0, 256},
        "LastPosX", "LastPosY");
    // Now for the experiment
    auto hTLExp =
        dfExpL1.Histo1D({"hExp", "TL distribution for exp;TL (mm);Counts", 200, 0, 200}, "MergerData.fLight.fTL");
    auto hLastVoxelExp = dfExpL1Positions.Histo2D(
        {"hLastVoxelExp", "Last voxel position for exp;LastPosX (mm);LastPosY (mm)", 256 / 2, 0, 256, 256 / 2, 0, 256},
        "LastPosX", "LastPosY");

    auto c1 = new TCanvas("c1", "c1", 1200, 600);
    c1->Divide(2, 2);
    c1->cd(1);
    hTLSimu->DrawClone();
    c1->cd(2);
    hTLExp->DrawClone();
    c1->cd(3);
    hLastVoxelSimu->DrawClone();
    c1->cd(4);
    hLastVoxelExp->DrawClone();

}