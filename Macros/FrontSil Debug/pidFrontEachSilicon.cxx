#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TH2D.h"
#include "TString.h"

#include <string>
#include <vector>

#include "../../PrettyStyle.C"


void pidFrontEachSilicon()
{
    PrettyStyle(true, true);

    std::string beam {"11Li"};
    std::string target {"d"};
    std::string light {"d"};

    // --------------------------------------------------
    // Get the data for all runs to see pid for each silicon
    // --------------------------------------------------

    std::string dataconf {};
    if(beam == "11Li")
        dataconf = "../../configs/data_11Li.conf";
    else if(beam == "7Li")
        dataconf = "../../configs/data_7Li.conf";
    else
        throw std::runtime_error("Beam cannot differ from 11Li or 7Li");

    // Read data
    ActRoot::DataManager dataman {dataconf, ActRoot::ModeType::EMerge};
    auto chain {dataman.GetChain()};

    // RDataFrame
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame dforigin {*chain};

    // Filter silicon pads
    auto df = dforigin.Filter(
        [](ActRoot::MergerData& m)
        {
            // Mask L0_9
            if(m.fLight.fNs.size())
                if((m.fLight.fLayers.front() == "l0"))
                    return false;
            // Mask F0_2
            if(m.fLight.fNs.size())
                if((m.fLight.fLayers.front() == "f0") && (m.fLight.fNs.front() == 2))
                    return false;
            // Mask R0_3
            if(m.fLight.fNs.size())
                if((m.fLight.fLayers.front() == "r0"))
                    return false;

            if(m.fRun > 35 && m.fRun < 45)
            {
                if(!m.fLight.fLayers.empty() && m.fLight.fLayers.front() == "f0")
                {
                    if(!m.fLight.fNs.empty() && m.fLight.fNs.front() == 5)
                    {
                        return false;
                    }
                    else
                        return true;
                }
                else
                    return true;
            }
            else
                return true;
        },
        {"MergerData"});

    // Only f0 events
    auto def = df.Filter([](ActRoot::MergerData& m)
                         { return m.fLight.GetNLayers() == 1 && m.fLight.GetLayer(0) == "f0"; }, {"MergerData"});


    // --------------------------------------------------
    // Create histograms (one per silicon)
    // --------------------------------------------------
    std::vector<TH2D*> hSil(12, nullptr);

    for(int i = 0; i < 12; i++)
    {
        hSil[i] = new TH2D(TString::Format("hPID_sil_%d", i),
                           TString::Format("PID silicon %d;E_{Sil} [MeV];Q_{ave}", i), 450, 0, 55, 60, 0, 500);
    }


    // --------------------------------------------------
    // Fill histograms
    // --------------------------------------------------
    def.Foreach(
        [&hSil](const ActRoot::MergerData& mer)
        {
            // Safety checks
            if(mer.fLight.fNs.empty())
                return;
            if(mer.fLight.fEs.empty())
                return;

            int sil = mer.fLight.fNs.front();
            if(sil < 0 || sil >= 12)
                return;

            hSil[sil]->Fill(mer.fLight.fEs.front(), mer.fLight.fQave);
        },
        {"MergerData"});


    // --------------------------------------------------
    // Draw
    // --------------------------------------------------
    auto* c = new TCanvas("c_pid_f0", "PID for f0", 1200, 800);
    c->Divide(4, 3);

    for(int i = 0; i < 12; i++)
    {
        c->cd(i + 1);
        hSil[i]->Draw("colz");
    }
}
