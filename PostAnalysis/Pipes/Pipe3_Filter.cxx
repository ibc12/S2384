


#ifndef PIPE3_FILTER_H
#define PIPE3_FILTER_H
#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActTPCData.h"

#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RResultPtr.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TColor.h"
#include "THStack.h"
#include "TString.h"
#include "TStyle.h"
#include "TVirtualPad.h"

#include <string>
#include <vector>

#include "../HistConfig.h"

// Filter the events reacting with CF4 that cannot be filtered with actroot -f
void Pipe3_Filter(const std::string& beam, const std::string& target, const std::string& light)
{
    std::string dataconf {};
    if(beam == "11Li")
        dataconf = "./../configs/data.conf";
    else if(beam == "7Li")
        dataconf = "./../configs/data_7Li.conf";
    else
        throw std::runtime_error("Beam cannot differ from 11Li or 7Li");

    // Read data
    ActRoot::DataManager dataman {dataconf, ActRoot::ModeType::EReadTPC};

    auto infile {TString::Format("./Outputs/tree_ex_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str())};

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"Final_Tree", infile.Data()};
    // Apply filter in average charge and reaches end of ACTAR
    auto dfFilter = df.Filter(
        [](ActRoot::MergerData& m)
        {
            if(m.fHeavy.fQave > 2000.) // Save value, max charge from 11Li/7Li is around 1800
                return false;
            return true;
        },
        {"MergerData"}).Filter(
            [](ActRoot::MergerData& m)
            {
                auto TL {m.fHeavy.fTL};
                auto theta {m.fThetaHeavy * TMath::DegToRad()};
                auto z_end {m.fRP.X() + TL * TMath::Cos(theta)};
                if(z_end < 240) // ACTAR ends at 256
                    return false;
                return true;
            },
            {"MergerData"});
    // Save!
    auto outfile {TString::Format("./Outputs/tree_ex_%s_%s_%s_filtered.root", beam.c_str(), target.c_str(), light.c_str())};
    dfFilter.Snapshot("Final_Tree", outfile);
    std::cout << "Saving Final_Tree in " << outfile << '\n';
}

#endif
