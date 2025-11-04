#include "ActMergerData.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TH1D.h"
#include "TString.h"

#include <fstream>
#include <iostream>
#include <string>

void DebugExElastic7Li()
{
    std::string beam {"11Li"};
    std::string target {"d"};
    std::string light {"p"};

    // Get file output from pipe2
    TString filename {TString::Format("../PostAnalysis/Outputs/tree_ex_%s_%s_%s_filtered.root", beam.c_str(),
                                      target.c_str(), light.c_str())};
    ROOT::RDataFrame df {"Final_Tree", filename};

    auto dfFiltered {df.Filter(
                           [](double& ex)
                           {
                               if(ex < -1)
                                   return true;
                               else
                                   return false;
                           },
                           {"Ex"})
                         .Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == true; }, {"MergerData"})};

    std::ofstream outFile("./Outputs/ExBellowNegative1_11Li_dp.dat");
    dfFiltered.Foreach([&](ActRoot::MergerData& m) { m.Stream(outFile); }, {"MergerData"});
    outFile.close();

    auto dfQMax = dfFiltered.Define("Qmax",
                                    [](const ActRoot::MergerData& m, const ActRoot::TPCData& tpc)
                                    {
                                        int heavyIdx = m.fHeavyIdx;
                                        float qmax = -1;
                                        for(auto& c : tpc.fClusters)
                                        {
                                            auto voxeles = c.GetVoxels();

                                            for(auto& v : voxeles)
                                            {
                                                if(v.GetCharge() > qmax)
                                                    qmax = v.GetCharge();
                                            }
                                        }
                                        return qmax;
                                    },
                                    {"MergerData", "GETTree_TPCData"});

    // Print the calculation of Qmax for debugging
    dfQMax.Foreach(
        [](ActRoot::MergerData& m, ActRoot::TPCData& tpc, float qmax)
        {
            std::cout << "--------------------------------" << std::endl;
            std::cout << "--------------------------------" << std::endl;
            std::cout << "--------------------------------" << std::endl;
            std::cout << "--------------------------------" << std::endl;
            int heavyIdx = m.fHeavyIdx;
            float computedQmax = -1;
            for(auto& c : tpc.fClusters)
            {
                auto voxeles = c.GetRefToVoxels();

                for(auto& v : voxeles)
                {
                    //std::cout << "Voxel Charge: " << v.GetCharge() << std::endl;
                    if(v.GetCharge() > computedQmax)
                        computedQmax = v.GetCharge();
                }
            }
            std::cout << "Qmax from Define: " << qmax << std::endl;
            std::cout << "Computed Qmax: " << computedQmax << std::endl;
            std::cout << "--------------------------------" << std::endl;
            std::cout << "--------------------------------" << std::endl;
            std::cout << "--------------------------------" << std::endl;
        },
        {"MergerData", "GETTree_TPCData", "Qmax"});

    // Plot Qave of fHeavy and Qmax and Ex
    auto hQaveHeavy = dfQMax.Histo1D({"hQaveHeavy", "Qave Heavy; Qave Heavy (a.u.); Counts", 100, 0, 5000},
                                     "MergerData.fHeavy.fQave");
    auto hQmax = dfQMax.Histo1D({"hQmax", "Qmax Heavy; Qmax Heavy (a.u.); Counts", 100, -10, 5000}, "Qmax");
    auto hEx = dfQMax.Histo1D({"hEx", "Ex; Ex (MeV); Counts", 100, -5, 10}, "Ex");
    // Qmax vs Ex
    auto hQmaxEx =
        dfQMax.Histo2D({"hQmaxEx", "Qmax vs Ex; Qmax (a.u.); Ex (MeV)", 100, -10, 5000, 100, -5, 10}, "Qmax", "Ex");
    TCanvas* c1 = new TCanvas("c1", "Qave Heavy", 800, 600);
    c1->DivideSquare(4);
    c1->cd(1);
    hQmax->DrawClone();
    c1->cd(2);
    hQaveHeavy->DrawClone();
    c1->cd(3);
    hEx->DrawClone();
    c1->cd(4);
    hQmaxEx->DrawClone();
}