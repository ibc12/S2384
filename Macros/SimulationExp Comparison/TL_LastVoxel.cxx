#include "ActInputParser.h"
#include "ActMergerData.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TString.h"

#include <string>
#include <vector>

constexpr int yMinExclusionZone = 55;
constexpr int yMaxExclusionZone = 70;

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
    TString fileName {TString::Format(
        "../../Simulation/Outputs/%s/test_charge_threshold/%s_%s_TRIUMF_Eex_%.3f_nPS_%d_pPS_%d_L1_2e6Thresh.root",
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
    auto dfExpL1Positions =
        dfExpL1
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
                        return lastVoxelPos.X() ; // Return in mm
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
                        return lastVoxelPos.Y(); // Return in mm
                    },
                    {"MergerData", "TPCData"})
            .Define("nVoxelsOutExclusionZone",
                    [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                    {
                        auto lightIdx = m.fLightIdx;
                        auto lightCl = tpc.fClusters[lightIdx];
                        auto voxels = lightCl.GetVoxels();

                        std::set<std::pair<int, int>> uniquePads;

                        for(const auto& voxel : voxels)
                        {
                            auto pos = voxel.GetPosition();
                            if(pos.Y() < yMinExclusionZone || pos.Y() > yMaxExclusionZone)
                                uniquePads.insert({static_cast<int>(pos.X()), static_cast<int>(pos.Y())});
                        }

                        return static_cast<int>(uniquePads.size());
                    },
                    {"MergerData", "TPCData"});

    // Define for the simulation the phi between -180 and 180 (now is 0 and 2pi)
    auto dfSimu = df.Define("phi3CMdeg",
                            [](double phi)
                            {
                                return (phi - M_PI) * 180.0 / M_PI; // Convert to degrees
                            },
                            {"phi3CM"});

    // Get the TL distribution for simu
    auto hTLSimu = dfSimu.Histo1D({"hSimu", "TL distribution for simu;TL (mm);Counts", 200, 0, 200}, "TL");
    auto hLastVoxelSimu = dfSimu.Histo2D({"hLastVoxelSimu", "Last voxel position for simu;LastPosX (mm);LastPosY (mm)",
                                          256 / 2, 0, 256, 256 / 2, 0, 256},
                                         "LastPosX", "LastPosY");
    // Now for the experiment
    auto hTLExp =
        dfExpL1.Histo1D({"hExp", "TL distribution for exp;TL (mm);Counts", 200, 0, 200}, "MergerData.fLight.fTL");
    auto hLastVoxelExp = dfExpL1Positions.Histo2D(
        {"hLastVoxelExp", "Last voxel position for exp;LastPosX (mm);LastPosY (mm)", 256 / 2, 0, 256 / 2, 256 / 2, 0, 256 / 2},
        "LastPosX", "LastPosY");

    // Get theta distributions for the experiment and simulation
    auto hThetaSimu =
        dfSimu.Histo1D({"hThetaSimu", "Theta distribution for simu;Theta (deg);Counts", 120, 0, 120}, "theta3Lab");
    auto hThetaExp = dfExpL1.Histo1D({"hThetaExp", "Theta distribution for exp;Theta (deg);Counts", 120, 0, 120},
                                     "MergerData.fThetaLight");

    // Get the phi distributions for the experiment and simulation
    auto hPhiSimu =
        dfSimu.Histo1D({"hPhiSimu", "Phi distribution for simu;Phi (deg);Counts", 360, -180, 180}, "phi3CMdeg");
    auto hPhiExp = dfExpL1.Histo1D({"hPhiExp", "Phi distribution for exp;Phi (deg);Counts", 360, -180, 180},
                                   "MergerData.fPhiLight");

    // Get the number of pads out of the exclusion zone for simu and exp
    auto hPadsSimu = dfSimu.Histo1D(
        {"hPadsSimu", "Number of pads out of exclusion zone for simu;N pads;Counts", 100, 0, 100}, "nPads");
    auto hPadsExp = dfExpL1Positions.Histo1D(
        {"hPadsExp", "Number of pads out of exclusion zone for exp;N pads;Counts", 100, 0, 100},
        "nVoxelsOutExclusionZone");

    // Get the number of pads as a function of phi for simu and exp
    auto hPadsPhiSimu =
        dfSimu.Histo2D({"hPadsPhiSimu", "Number of pads out of exclusion zone vs phi for simu;Phi (deg);N pads", 360,
                        -180, 180, 100, 0, 100},
                       "phi3CMdeg", "nPads");
    auto hPadsPhiExp = dfExpL1Positions.Histo2D(
        {"hPadsPhiExp", "Number of pads out of exclusion zone vs phi for exp;Phi (deg);;N pads", 360, -180, 180, 100, 0,
         100},
        "MergerData.fPhiLight", "nVoxelsOutExclusionZone");

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

    auto c2 = new TCanvas("c2", "c2", 1200, 600);
    c2->Divide(2, 2);
    c2->cd(1);
    hThetaSimu->DrawClone();
    c2->cd(2);
    hThetaExp->DrawClone();
    c2->cd(3);
    hPhiSimu->DrawClone();
    c2->cd(4);
    hPhiExp->DrawClone();

    auto c3 = new TCanvas("c3", "c3", 1200, 600);
    c3->Divide(2, 2);
    c3->cd(1);
    hPadsSimu->DrawClone();
    c3->cd(2);
    hPadsExp->DrawClone();
    c3->cd(3);
    hPadsPhiSimu->DrawClone("colz");
    c3->cd(4);
    hPadsPhiExp->DrawClone("colz");
}