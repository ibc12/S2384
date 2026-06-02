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
    std::string beam {"7Li"};
    std::string target {"d"};
    std::string light {"d"};

    // Get file output from pipe2
    TString filename {TString::Format("../PostAnalysis/Outputs/tree_ex_F_%s_%s_%s_filtered.root", beam.c_str(),
                                      target.c_str(), light.c_str())};
    ROOT::RDataFrame df {"Final_Tree", filename};

    auto dfFiltered {df.Filter(
                           [](double& ex)
                           {
                               if(ex > 1 && ex < 2)
                                   return true;
                               else
                                   return false;
                           },
                           {"Ex"})
                         .Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == true; }, {"MergerData"})};

    std::ofstream outFile("./Outputs/ExBetween1And2_7Li_dd.dat");
    dfFiltered.Foreach([&](ActRoot::MergerData& m) { m.Stream(outFile); }, {"MergerData"});
    outFile.close();
}