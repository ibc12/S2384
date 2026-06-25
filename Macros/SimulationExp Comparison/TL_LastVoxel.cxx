#include "ActInputParser.h"
#include "ActMergerData.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TROOT.h"
#include "TString.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

constexpr int yMinExclusionZone = 56;
constexpr int yMaxExclusionZone = 71;

void TL_LastVoxel()
{
    // Get data for 7Li (d,d) elastic scattering
    std::string beam {"11Li"};
    std::string target {"2H"};
    std::string targetExp {"d"};
    std::string light {"2H"};
    std::string lightExp {"d"};
    std::string heavy {"7Li"};
    double Ex {0.}; // MeV
    // Get simu data
    TString fileName {TString::Format(
        "../../Simulation/Outputs/%s/test_charge_threshold/%s_%s_TRIUMF_Eex_%.3f_nPS_%d_pPS_%d_L1_1e6Thresh.root",
        beam.c_str(), target.c_str(), light.c_str(), Ex, 0, 0)};
    // Get exp data
    TString fileNameExp {TString::Format("../../PostAnalysis/Outputs/tree_ex_F_%s_%s_%s.root", beam.c_str(),
                                         targetExp.c_str(), lightExp.c_str())};

    ROOT::EnableImplicitMT();

    auto df = ROOT::RDataFrame("SimulationTTree", fileName);
    auto dfExp = ROOT::RDataFrame("Final_Tree", fileNameExp);

    // Filter exp data to get only L1 events
    auto dfExpL1 = dfExp.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"})
                       .Filter([](ActRoot::MergerData& m) { return m.fRun ; }, {"MergerData"});
    // auto dfExpL1 = dfExpL1.Filter([](double & ex) { return ex < 1; }, {"Ex"});

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
                        return lastVoxelPos.X(); // Return in pad
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
                        return lastVoxelPos.Y(); // Return in pad
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
                    {"MergerData", "TPCData"})
            .Define("nVoxelsRawOutExclusionZone",
                    [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                    {
                        auto lightIdx = m.fLightIdx;
                        auto lightCl = tpc.fClusters[lightIdx];
                        auto voxels = tpc.fRaw;

                        std::set<std::pair<int, int>> uniquePads;

                        for(const auto& voxel : voxels)
                        {
                            auto pos = voxel.GetPosition();
                            if(pos.Y() < yMinExclusionZone || pos.Y() > yMaxExclusionZone)
                                uniquePads.insert({static_cast<int>(pos.X()), static_cast<int>(pos.Y())});
                        }

                        return static_cast<int>(uniquePads.size());
                    },
                    {"MergerData", "TPCData"})
            .Define("nVoxelsBeamOutExclusionZone",
                    [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                    {
                        auto beamIdx = m.fBeamIdx;
                        auto beamCl = tpc.fClusters[beamIdx];
                        auto voxels = beamCl.GetVoxels();

                        std::set<std::pair<int, int>> uniquePads;

                        for(const auto& voxel : voxels)
                        {
                            auto pos = voxel.GetPosition();
                            if(pos.Y() < yMinExclusionZone || pos.Y() > yMaxExclusionZone)
                                uniquePads.insert({static_cast<int>(pos.X()), static_cast<int>(pos.Y())});
                        }

                        return static_cast<int>(uniquePads.size());
                    },
                    {"MergerData", "TPCData"})
            .Define("nVoxelsHeavyOutExclusionZone",
                    [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                    {
                        auto heavyIdx = m.fHeavyIdx;
                        auto heavyCl = tpc.fClusters[heavyIdx];
                        auto voxels = heavyCl.GetVoxels();

                        std::set<std::pair<int, int>> uniquePads;

                        for(const auto& voxel : voxels)
                        {
                            auto pos = voxel.GetPosition();
                            if(pos.Y() < yMinExclusionZone || pos.Y() > yMaxExclusionZone)
                                uniquePads.insert({static_cast<int>(pos.X()), static_cast<int>(pos.Y())});
                        }

                        return static_cast<int>(uniquePads.size());
                    },
                    {"MergerData", "TPCData"})
            .Define("nVoxelsTotal",
                    [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                    {
                        auto lightIdx = m.fLightIdx;
                        auto heavyIdx = m.fHeavyIdx;
                        auto beamIdx = m.fBeamIdx;
                        auto lightCl = tpc.fClusters[lightIdx];
                        auto heavyCl = tpc.fClusters[heavyIdx];
                        auto beamCl = tpc.fClusters[beamIdx];
                        auto voxelsLight = lightCl.GetVoxels();
                        auto voxelsHeavy = heavyCl.GetVoxels();
                        auto voxelsBeam = beamCl.GetVoxels();

                        std::set<std::pair<int, int>> uniquePads;

                        for(const auto& voxel : voxelsLight)
                        {
                            auto pos = voxel.GetPosition();
                            if(pos.Y() < yMinExclusionZone || pos.Y() > yMaxExclusionZone)
                                uniquePads.insert({static_cast<int>(pos.X()), static_cast<int>(pos.Y())});
                        }
                        for(const auto& voxel : voxelsHeavy)
                        {
                            auto pos = voxel.GetPosition();
                            if(pos.Y() < yMinExclusionZone || pos.Y() > yMaxExclusionZone)
                                uniquePads.insert({static_cast<int>(pos.X()), static_cast<int>(pos.Y())});
                        }
                        for(const auto& voxel : voxelsBeam)
                        {
                            auto pos = voxel.GetPosition();
                            if(pos.Y() < yMinExclusionZone || pos.Y() > yMaxExclusionZone)
                                uniquePads.insert({static_cast<int>(pos.X()), static_cast<int>(pos.Y())});
                        }
                        for(const auto& voxel : tpc.fRaw)
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
    auto hLastVoxelExp =
        dfExpL1Positions.Histo2D({"hLastVoxelExp", "Last voxel position for exp;LastPosX (mm);LastPosY (mm)", 256 / 2,
                                  0, 256 / 2, 256 / 2, 0, 256 / 2},
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
        {"hPadsSimu", "Number of pads out of exclusion zone for simu;N pads;Counts", 150, 0, 150}, "nPads");
    auto hPadsExp = dfExpL1Positions.Histo1D(
        {"hPadsExp", "Number of pads out of exclusion zone for exp;N pads;Counts", 150, 0, 150},
        "nVoxelsOutExclusionZone");

    // Get the number of pads as a function of phi for simu and exp
    auto hPadsPhiSimu =
        dfSimu.Histo2D({"hPadsPhiSimu", "Number of pads out of exclusion zone vs phi for simu;Phi (deg);N pads", 360,
                        -180, 180, 150, 0, 150},
                       "phi3CMdeg", "nPads");
    auto hPadsPhiExp = dfExpL1Positions.Histo2D(
        {"hPadsPhiExp", "Number of pads out of exclusion zone vs phi for exp;Phi (deg);;N pads", 360, -180, 180, 150, 0,
         150},
        "MergerData.fPhiLight", "nVoxelsOutExclusionZone");

    // Get the number of pads as a function of the TL for simu and exp
    auto hPadsTLSimu = dfSimu.Histo2D(
        {"hPadsTLSimu", "Number of pads out of exclusion zone vs TL for simu;TL (mm);N pads", 200, 0, 200, 150, 0, 150},
        "TL", "nPads");
    auto hPadsTLExp = dfExpL1Positions.Histo2D(
        {"hPadsTLExp", "Number of pads out of exclusion zone vs TL for exp;TL (mm);N pads", 200, 0, 200, 150, 0, 150},
        "MergerData.fLight.fTL", "nVoxelsOutExclusionZone");

    // Get the numebr of pads for the raw voxels (not clustered) for the exp data
    auto hPadsRawExpPhi = dfExpL1Positions.Histo2D(
        {"hPadsRawExpPhi", "Number of pads out of exclusion zone for raw voxels for exp;Phi (deg);N pads", 360, -180,
         180, 150, 0, 150},
        "MergerData.fPhiLight", "nVoxelsRawOutExclusionZone");
    auto hPadsRawExpTL = dfExpL1Positions.Histo2D(
        {"hPadsRawExpTL", "Number of pads out of exclusion zone for raw voxels for exp;TL (mm);N pads", 200, 0, 200,
         150, 0, 150},
        "MergerData.fLight.fTL", "nVoxelsRawOutExclusionZone");
    // Now the heavy cluster
    auto hPadsHeavyExpPhi = dfExpL1Positions.Histo2D(
        {"hPadsHeavyExpPhi", "Number of pads out of exclusion zone for heavy cluster for exp;Phi (deg);N pads", 360,
         -180, 180, 150, 0, 150},
        "MergerData.fPhiLight", "nVoxelsHeavyOutExclusionZone");
    auto hPadsHeavyExpTL = dfExpL1Positions.Histo2D(
        {"hPadsHeavyExpTL", "Number of pads out of exclusion zone for heavy cluster for exp;TL (mm);N pads", 200, 0,
         200, 150, 0, 150},
        "MergerData.fLight.fTL", "nVoxelsHeavyOutExclusionZone");
    // Now the beam cluster
    auto hPadsBeamExpPhi = dfExpL1Positions.Histo2D(
        {"hPadsBeamExpPhi", "Number of pads out of exclusion zone for beam cluster for exp;Phi (deg);N pads", 360, -180,
         180, 150, 0, 150},
        "MergerData.fPhiLight", "nVoxelsBeamOutExclusionZone");
    auto hPadsBeamExpTL = dfExpL1Positions.Histo2D(
        {"hPadsBeamExpTL", "Number of pads out of exclusion zone for beam cluster for exp;TL (mm);N pads", 200, 0, 200,
         150, 0, 150},
        "MergerData.fLight.fTL", "nVoxelsBeamOutExclusionZone");

    // Get the plots pads vs TL and Phi against nPads total out of the exclusion zone (light + heavy + beam + raw
    // voxels)
    auto hPadsAllLightExpPhi = dfExpL1Positions.Histo2D(
        {"hPadsAllLightExpPhi", "Number of pads out of exclusion zone for all + light cluster for exp;Phi (deg);N pads",
         360, -180, 180, 150, 0, 150},
        "MergerData.fPhiLight", "nVoxelsTotal");
    auto hPadsAllLightExpTL = dfExpL1Positions.Histo2D(
        {"hPadsAllLightExpTL", "Number of pads out of exclusion zone for all + light cluster for exp;TL (mm);N pads",
         200, 0, 200, 150, 0, 150},
        "MergerData.fLight.fTL", "nVoxelsTotal");
    auto hPadsAll = dfExpL1Positions.Histo1D(
        {"hPadsAll", "Number of pads out of exclusion zone for all + light cluster for exp;N pads;Counts", 150, 0, 150},
        "nVoxelsTotal");

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

    auto c4 = new TCanvas("c4", "c4", 1200, 600);
    c4->Divide(2, 2);
    c4->cd(1);
    hPadsTLSimu->DrawClone("colz");
    c4->cd(2);
    hPadsTLExp->DrawClone("colz");
    c4->cd(3);
    hPadsRawExpPhi->DrawClone("colz");
    c4->cd(4);
    hPadsRawExpTL->DrawClone("colz");

    auto c5 = new TCanvas("c5", "c5", 1200, 600);
    c5->Divide(2, 2);
    c5->cd(1);
    hPadsHeavyExpPhi->DrawClone("colz");
    c5->cd(2);
    hPadsBeamExpPhi->DrawClone("colz");
    c5->cd(3);
    hPadsRawExpPhi->DrawClone("colz");
    c5->cd(4);
    hPadsPhiExp->DrawClone("colz");

    auto c6 = new TCanvas("c6", "c6", 1200, 600);
    c6->Divide(2, 2);
    c6->cd(1);
    hPadsHeavyExpTL->DrawClone("colz");
    c6->cd(2);
    hPadsBeamExpTL->DrawClone("colz");
    c6->cd(3);
    hPadsRawExpTL->DrawClone("colz");
    c6->cd(4);
    hPadsTLExp->DrawClone("colz");

    auto c7 = new TCanvas("c7", "c7", 1200, 600);
    c7->Divide(2, 2);
    c7->cd(1);
    hPadsAllLightExpPhi->DrawClone("colz");
    c7->cd(2);
    hPadsAllLightExpTL->DrawClone("colz");
    c7->cd(3);
    hPadsAll->DrawClone();


    // Get some events to inspect
    // auto dfGetEvents = dfExpL1Positions.Filter([](int nPads) { return nPads > 70; }, {"nVoxelsOutExclusionZone"});
    // std::ofstream outEvent("./Outputs/nPads100.dat");
    // dfGetEvents.Foreach(
    //     [&outEvent](ActRoot::MergerData& m)
    //     {
    //         m.Stream(outEvent);
    //     },
    //     {"MergerData"});
}