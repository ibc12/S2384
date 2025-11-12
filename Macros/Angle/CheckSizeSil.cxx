#include "ActMergerData.h"
#include "ActSilMatrix.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TString.h"

#include <filesystem>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

#include "../SilVetos/GetContourFuncs.cxx" // For FindBestFit, FitToContour, etc.

using namespace ROOT;

std::pair<double, double> FitAngularSize(TH1D* h)
{
    double w {5.0};
    double s {0.5};
    auto* func {FindBestFit(h, w, s)};
    h = ScaleWithFunc(h, func);
    double thresh {0.65};
    double width {10.0};
    return FitToCountour(h, thresh, width);
}

void CheckSizeSil()
{
    // ROOT::EnableImplicitMT();
    // Read data
    std::string beam {"7Li"};
    std::string target {"d"};
    std::string light {"p"};
    auto filename {TString::Format("../../PostAnalysis/Outputs/tree_ex_%s_%s_%s.root", beam.c_str(), target.c_str(),
                                   light.c_str())};
    ROOT::RDataFrame df {"Final_Tree", filename};

    std::string layer {"l0"};
    std::vector<int> indexes;
    if(layer == "l0")
        indexes = {0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11};
    // indexes = {5, 8};
    else if(layer == "f0")
        indexes = {4, 7};
    else if(layer == "r0")
        indexes = {0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11};

    double distance {};
    if(layer == "r0")
        distance = 71.; // mm
    else if(layer == "l0")
        distance = 59.; // mm
    distance += 128.; // From the middle of ACTAR to end of pad plane

    // Prepare output folder
    std::filesystem::create_directories("./Outputs/Angles");

    // Histograms by silicon ID
    std::map<int, std::shared_ptr<TH1D>> hAngles;

    // Create one histogram per silicon index (1–8, adjust if needed)
    for(int idx : indexes)
    {
        TString hname = TString::Format("hTheta_sil%d", idx);
        TString htitle = TString::Format("Angular distribution for silicon %d", idx);
        hAngles[idx] = std::make_shared<TH1D>(hname, htitle, 180, 0, 180);
    }

    // Fill histograms from dataframe
    df.Foreach(
        [&](const ActRoot::MergerData& m)
        {
            if(m.fLight.fLayers.empty() || m.fLight.fNs.empty())
                return;

            const auto& lyr = m.fLight.fLayers.front();
            const auto& n = m.fLight.fNs.front();

            if(lyr == layer)
            {
                auto it = hAngles.find(n);
                if(it != hAngles.end() && it->second) // ✅ check existence and non-null
                    it->second->Fill(m.fThetaLight);
            }
        },
        {"MergerData"});

    // 1. Print number of entries per silicon
    std::cout << "\n=== Histogram summary ===\n";
    for(auto& [idx, h] : hAngles)
        std::cout << "Silicon " << idx << ": " << h->GetEntries() << " entries\n";

    for(auto& [idx, h] : hAngles)
    {
        auto [amin, amax] = FitAngularSize(h.get());
        std::cout << "Silicon " << idx << ": min = " << amin << ", max = " << amax << " deg\n";

        // Convert degrees to radians
        double amin_rad = amin * M_PI / 180.0;
        double amax_rad = amax * M_PI / 180.0;

        // Physical size
        double size_mm = distance * (tan(amax_rad) - tan(amin_rad));
        std::cout << "Silicon " << idx << ": estimated size ~ " << size_mm << " mm\n";
    }

    // Draw histograms
    auto* c1 = new TCanvas("c1", "Angular distributions", 1200, 800);
    c1->DivideSquare(hAngles.size());
    int pad = 1;
    for(auto& [idx, h] : hAngles)
    {
        c1->cd(pad++);
        h->DrawClone();
    }


    // Save all histograms
    TString outFile = TString::Format("./Outputs/histos_all_silicons_%s.root", layer.c_str());
    auto fOut = std::make_unique<TFile>(outFile, "RECREATE");
    for(auto& [idx, h] : hAngles)
        h->Write();
    fOut->Close();
}
