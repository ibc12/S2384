#include "ActMergerData.h"

#include "ROOT/RDataFrame.hxx"

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
    TString filename {
        TString::Format("../PostAnalysis/Outputs/tree_ex_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str())};
    ROOT::RDataFrame df {"Final_Tree", filename};

    auto dfFiltered {df.Filter(
                           [](double& ex)
                           {
                               if(ex < -1.0)
                                   return true;
                               else
                                   return false;
                           },
                           {"Ex"})
                         .Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == true; }, {"MergerData"})};

    std::ofstream outFile("./Outputs/ExBellowNegative1_7Li_Elastic.dat");
    dfFiltered.Foreach([&](ActRoot::MergerData& m) { m.Stream(outFile); }, {"MergerData"});
    outFile.close();
}