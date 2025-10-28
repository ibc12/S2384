#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"

#include <string>
#include <vector>
#include <iostream>

void ExSpectra()
{
    // Read output from simulation and plot all the peaks you want
    std::string beam {"7Li"};
    std::string target {"2H"};
    std::string light {"2H"};
    std::vector<double> Exs {0.0, 0.5};
    int neutronPS {0};
    int protonPS {0};


    auto c = new TCanvas("canvas", "Eex Spectrum", 800, 600);
    auto hEx {new TH1D("hEx", "Excitation Energy Spectrum;E_{ex} (MeV);Counts", 100, -1, 3)};
    // Open files and plot Eex branch from TTree SimulationTTree
    for(const auto& Ex : Exs)
    {
        TString fileName {TString::Format("../Outputs/%s/%s_%s_TRIUMF_Eex_%.3f_nPS_%d_pPS_%d.root", beam.c_str(),
                                          target.c_str(), light.c_str(), Ex, neutronPS, protonPS)};
        auto* inFile {TFile::Open(fileName, "read")};
        if(!inFile || inFile->IsZombie())
        {
            std::cout << "Error opening file: " << fileName.Data() << std::endl;
            continue;
        }
        auto* inTree {static_cast<TTree*>(inFile->Get("SimulationTTree"))};
        if(!inTree)
        {
            std::cout << "Error: TTree 'SimulationTTree' not found in file: " << fileName.Data() << std::endl;
            inFile->Close();
            continue;
        }
        // Plot branch Eex adding to total hEx
        double Eex_value {};
        inTree->SetBranchAddress("Eex", &Eex_value);
        int nEntries = inTree->GetEntries();
        for(int i = 0; i < nEntries; ++i)
        {
            inTree->GetEntry(i);
            hEx->Fill(Eex_value);
        }
        inFile->Close();
    }
    hEx->Draw();
}