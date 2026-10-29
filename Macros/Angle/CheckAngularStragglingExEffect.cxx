#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TTree.h"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

void CheckAngularStragglingExEffect()
{
    // -------------------------------------------------------------------------
    // Configuration
    // -------------------------------------------------------------------------
    const std::string basePath = "../../Simulation/Outputs/7Li/test_ang_straggling_L1/";
    const std::vector<std::string> angStrValues = {"3-5", "4", "4-5", "5", "5-5", "6", "6-5"};
    const int nBins = 100;
    const double eMin = -1.0;
    const double eMax = 1.0;
    const bool fitGaussian = false;
    // -------------------------------------------------------------------------
    // ROOT style
    // -------------------------------------------------------------------------
    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(1111);
    // -------------------------------------------------------------------------
    // Canvas
    // -------------------------------------------------------------------------
    auto* canvas = new TCanvas("cAngularStraggling", "Angular straggling effect on excitation energy", 1400, 1000);
    canvas->Divide(2, 4, 0.01, 0.01);
    // Keep histograms alive after the function returns.
    std::vector<TH1D*> histograms;
    // -------------------------------------------------------------------------
    // Loop over angular straggling values
    // -------------------------------------------------------------------------
    for(std::size_t i = 0; i < angStrValues.size(); ++i)
    {
        const auto& angStr = angStrValues[i];

        const std::string fileName = basePath + "2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_" + angStr + "AngStr.root";

        std::cout << "\nProcessing angular straggling = " << angStr << '\n' << "File: " << fileName << '\n';

        // Open input file.
        std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "READ"));
        if(!file || file->IsZombie())
        {
            std::cerr << "ERROR: Cannot open file: " << fileName << '\n';
            continue;
        }
        // Retrieve simulation tree.
        auto* tree = dynamic_cast<TTree*>(file->Get("SimulationTTree"));
        if(!tree)
        {
            std::cerr << "ERROR: SimulationTTree not found in " << fileName << '\n';
            continue;
        }
        // Create a uniquely named histogram (it is attached to the open file for now,
        // which is what tree->Draw("Eex>>name") needs in order to find it).
        const std::string histName = "hEex_" + std::to_string(i);
        auto* hEex = new TH1D(histName.c_str(), ("Angular straggling = " + angStr + ";E_{ex} [MeV];Counts").c_str(),
                              nBins, eMin, eMax);
        hEex->SetLineColor(kBlue + 1);
        hEex->SetLineWidth(2);

        // Fill histogram from the simulation tree.
        const Long64_t nSelected = tree->Draw(("Eex>>" + histName).c_str(), "", "goff");
        if(nSelected < 0)
        {
            std::cerr << "ERROR: Failed to draw Eex from " << fileName << '\n';

            delete hEex;
            continue;
        }

        // Detach AFTER filling so the histogram survives the closing of the TFile.
        hEex->SetDirectory(nullptr);
        histograms.push_back(hEex);

        std::cout << "Entries in tree: " << tree->GetEntries() << "\nEntries selected: " << nSelected
                  << "\nEntries in histogram: " << hEex->GetEntries() << "\nMean: " << hEex->GetMean() << " MeV"
                  << "\nStandard deviation: " << hEex->GetStdDev() << " MeV\n";

        // Draw histogram.
        canvas->cd(i + 1);
        hEex->Draw("HIST");
        // Optional Gaussian fit.
        if(fitGaussian && hEex->GetEntries() > 0)
        {
            TF1* gaussian = new TF1(("fGaus_" + std::to_string(i)).c_str(), "gaus", eMin, eMax);
            const double mean = hEex->GetMean();
            const double sigma = hEex->GetStdDev();
            gaussian->SetParameters(hEex->GetMaximum(), mean, sigma > 0.0 ? sigma : 0.1);
            const int fitStatus = hEex->Fit(gaussian, "RQ");
            if(fitStatus == 0)
            {
                std::cout << "Gaussian fit:"
                          << "\n  Mean = " << gaussian->GetParameter(1) << " MeV"
                          << "\n  Sigma = " << gaussian->GetParameter(2) << " MeV"
                          << "\n  FWHM = " << 2.35482 * gaussian->GetParameter(2) << " MeV\n";
            }
        }
    }

    canvas->Modified();
    canvas->Update();

    std::cout << "\nComparison completed.\n";
}