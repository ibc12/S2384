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

#include "../PrettyStyle.C"

void checkDataML()
{
    std::string infile {"../PostAnalysis/Outputs/tree_preprocess_F_7Li.root"};
    std::string infile1 {"../PostAnalysis/Outputs/tree_preprocess_7Li.root"};

    // RDataFrame
    // ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"PreProcessed_Tree", infile};
    ROOT::RDataFrame df1 {"PreProcessed_Tree", infile1};

    // Get drift parameter from config file
    ActRoot::InputParser parser {};
    parser.ReadFile("../configs/detector.conf");
    auto driftBlock = parser.GetBlock("Merger");
    auto driftFactor = driftBlock->GetDouble("DriftFactor"); // in mm^2/us

    auto dfFiltered =
        df.Filter([](ActRoot::MergerData& mer, ActRoot::ModularData& mod)
                  { return (mod.Get("GATCONF") == 8) && (mer.fLightIdx != -1); }, {"MergerData", "ModularData"});
    auto df1Filtered =
        df1.Filter([](ActRoot::MergerData& mer, ActRoot::ModularData& mod)
                   { return (mod.Get("GATCONF") == 8) && (mer.fLightIdx != -1); }, {"MergerData", "ModularData"});

    auto dfDefines =
        dfFiltered
            .Define("X",
                    [](ActRoot::MergerData& mer, ActRoot::TPCData& tpc)
                    {
                        ROOT::RVec<float> x;
                        if(mer.fLightIdx < 0 || static_cast<std::size_t>(mer.fLightIdx) >= tpc.fClusters.size())
                            return x;
                        auto& cluster = tpc.fClusters.at(mer.fLightIdx);
                        auto& voxels = cluster.GetRefToVoxels();

                        for(const auto& v : voxels)
                        {
                            auto pos = v.GetPosition();
                            x.push_back(pos.X() * 2);
                        }
                        return x;
                    },
                    {"MergerData", "TPCData"})
            .Define("Y",
                    [](ActRoot::MergerData& mer, ActRoot::TPCData& tpc)
                    {
                        ROOT::RVec<float> y;
                        if(mer.fLightIdx < 0 || static_cast<std::size_t>(mer.fLightIdx) >= tpc.fClusters.size())
                            return y;
                        auto& cluster = tpc.fClusters.at(mer.fLightIdx);
                        auto& voxels = cluster.GetRefToVoxels();

                        for(const auto& v : voxels)
                        {
                            auto pos = v.GetPosition();
                            y.push_back(pos.Y() * 2);
                        }
                        return y;
                    },
                    {"MergerData", "TPCData"})
            .Define("Z",
                    [&driftFactor](ActRoot::MergerData& mer, ActRoot::TPCData& tpc)
                    {
                        ROOT::RVec<float> z;
                        if(mer.fLightIdx < 0 || static_cast<std::size_t>(mer.fLightIdx) >= tpc.fClusters.size())
                            return z;
                        auto& cluster = tpc.fClusters.at(mer.fLightIdx);
                        auto& voxels = cluster.GetRefToVoxels();

                        for(const auto& v : voxels)
                        {
                            auto pos = v.GetPosition();
                            z.push_back(pos.Z() * driftFactor);
                        }
                        return z;
                    },
                    {"MergerData", "TPCData"})
            .Define("Charge",
                    [](ActRoot::MergerData& mer, ActRoot::TPCData& tpc)
                    {
                        ROOT::RVec<float> charge;
                        if(mer.fLightIdx < 0 || static_cast<std::size_t>(mer.fLightIdx) >= tpc.fClusters.size())
                            return charge;
                        auto& cluster = tpc.fClusters.at(mer.fLightIdx);
                        auto& voxels = cluster.GetRefToVoxels();

                        for(const auto& v : voxels)
                            charge.push_back(v.GetCharge());
                        return charge;
                    },
                    {"MergerData", "TPCData"});


    std::string dataMLPath = "./Outputs/dataForML.root";
    // Get the entry 773 of the root file (x,y,z,charge structure of event)
    TFile* inFile = TFile::Open(dataMLPath.c_str(), "READ");
    if(!inFile || inFile->IsZombie())
    {
        std::cerr << "Error opening file: " << dataMLPath << std::endl;
        return;
    }
    TTree* tree = dynamic_cast<TTree*>(inFile->Get("ML_Tree"));
    if(!tree)
    {
        std::cerr << "Error: TTree 'ML_Tree' not found in file: " << dataMLPath << std::endl;
        inFile->Close();
        return;
    }
    std::vector<float>* x = nullptr;
    std::vector<float>* y = nullptr;
    std::vector<float>* z = nullptr;
    std::vector<float>* charge = nullptr;

    tree->SetBranchAddress("X", &x);
    tree->SetBranchAddress("Y", &y);
    tree->SetBranchAddress("Z", &z);
    tree->SetBranchAddress("Charge", &charge);
    tree->GetEntry(769);

    // Check what event of the df has the same charges as the vector charge
    auto dfEvent = dfDefines.Filter(
        [&](ActRoot::MergerData& mer, ActRoot::TPCData& tpc)
        {
            if(mer.fLightIdx < 0 || static_cast<std::size_t>(mer.fLightIdx) >= tpc.fClusters.size())
                return false;
            auto& cluster = tpc.fClusters.at(mer.fLightIdx);
            auto& voxels = cluster.GetRefToVoxels();
            if(voxels.size() != charge->size())
                return false;
            for(size_t i = 0; i < voxels.size(); ++i)
            {
                if(std::abs(voxels[i].GetCharge() - (*charge)[i]) > 1e-3)
                    return false;
            }
            return true;
        },
        {"MergerData", "TPCData"});
    auto count = dfEvent.Count();
    std::cout << "Number of events with same charges: " << *count << std::endl;
    dfEvent.Foreach(
        [](ActRoot::MergerData& mer)
        {
            std::cout << "Entry: " << mer.fEntry << ", Run: " << mer.fRun << std::endl;
        },
        {"MergerData"});

    // Plot the projections xy yz xz of the voxels for this event, with the charge as color
    auto* c = new TCanvas("cDataML", "Data for ML", 1200, 400);
    c->Divide(3, 1);
    c->cd(1);
    auto* hXY = new TH2D("hXY", "XY projection;X [mm];Y [mm]", 128, 0, 256, 128, 0, 256);
    hXY->SetOption("colz");
    for(size_t i = 0; i < x->size(); ++i)
    {
        hXY->Fill((*x)[i], (*y)[i], (*charge)[i]);
    }
    hXY->DrawClone("colz");
    c->cd(2);
    auto* hYZ = new TH2D("hYZ", "YZ projection;Y [mm];Z [mm]", 128, 0, 256, 200, 0, 400);
    hYZ->SetOption("colz");
    for(size_t i = 0; i < y->size(); ++i)
    {
        hYZ->Fill((*y)[i], (*z)[i], (*charge)[i]);
    }
    hYZ->DrawClone("colz");
    c->cd(3);
    auto* hXZ = new TH2D("hXZ", "XZ projection;X [mm];Z [mm]", 128, 0, 256, 200, 0, 400);
    hXZ->SetOption("colz");
    for(size_t i = 0; i < x->size(); ++i)
    {
        hXZ->Fill((*x)[i], (*z)[i], (*charge)[i]);
    }
    hXZ->DrawClone("colz");
}