#include "ActDataManager.h"
#include "ActSilData.h"
#include "ActTPCData.h"
#include "ActModularData.h"
#include "ActModularData.h"
#include "ActMergerData.h"
#include "ActTPCData.h"
#include "ActCluster.h"
#include "ActVoxel.h"

#include <ROOT/RDataFrame.hxx>

#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"
#include "TF1.h"

void gainVSmeshV()
{
    std::vector<TH1D*> chargeHistograms;
    std::vector<double> maxCharges {};
    for (int run = 124; run <= 128; ++run)
    {
        // build charge distribution histogram
        auto hist {new TH1D(("chargeDist_Run" + std::to_string(run)).c_str(),
                                    ("Charge Distribution Run " + std::to_string(run)).c_str(),
                                    100, 0, 5000)};
        // Construir nombre del archivo
        std::string filename = "../RootFiles/Cluster/Clusters_Run_0" + std::to_string(run) + ".root";

        // Abrir archivo
        TFile file(filename.c_str(), "READ");
        if (file.IsZombie()) {
            std::cerr << "No se pudo abrir " << filename << "\n";
            continue;
        }

        // Obtener TTree
        TTree *tree = nullptr;
        file.GetObject("GETTree", tree);
        if (!tree) {
            std::cerr << "No se encontró TTree en " << filename << "\n";
            continue;
        }

        // Conectar branch de TPCData
        ActRoot::TPCData *tpcData = nullptr;
        tree->SetBranchAddress("TPCData", &tpcData);

        double sumCharge = 0;
        Long64_t nEvents = tree->GetEntries();
        std::cout << "Procesando run " << run << " con " << nEvents << " eventos.\n";
        int charge0events = 0;

        // Recorrer eventos
        for (Long64_t i = 0; i < tree->GetEntries(); ++i)
        {
            tree->GetEntry(i);

            double maxCharge = 0;
            if (tpcData->fClusters.size() == 1)
            {
                auto voxels = tpcData->fClusters.front().GetRefToVoxels();
                for (const auto &v : voxels)
                {
                    if (v.GetIsSaturated())
                        continue;
                    if (v.GetCharge() > maxCharge)
                        maxCharge = v.GetCharge();
                }
            }
            if(maxCharge == 0)
            {
                charge0events++;
                nEvents = nEvents - 1;
            }
            else
                hist->Fill(maxCharge);
            sumCharge += maxCharge;
        }

        chargeHistograms.push_back(hist);

        double meanCharge = (nEvents > 0) ? sumCharge / nEvents : 0;
        maxCharges.push_back(meanCharge);
        std::cout << "Run " << run << " -> media maxCharge = " << meanCharge << " charge 0 events " << charge0events <<"\n";
    }

    std::vector<int> meshVoltages = {540, 530, 520, 510, 500};
    TGraph *graph = new TGraph(meshVoltages.size());
    for (size_t i = 0; i < meshVoltages.size(); ++i)
    {
        graph->SetPoint(i, meshVoltages[i], TMath::Log10(maxCharges[i]));
    }
    graph->SetTitle("Gain vs Mesh Voltage;Mesh Voltage [V];Log Mean Max Charge [a.u.]");
    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(kBlue);
    graph->SetMarkerSize(1.2);
    // Fit the data
    graph->Fit("pol1", "0Q"); // Fit with a polynomial of degree 2
    TF1 *fitFunc = graph->GetFunction("pol1");
    // Plot graph
    TCanvas *c1 = new TCanvas("c1", "Gain vs Mesh Voltage", 800, 600);
    graph->DrawClone("AP");
    fitFunc->SetLineColor(kRed);
    fitFunc->Draw("same");

    // Plot charge hists 
    TCanvas *c2 = new TCanvas("c2", "Charge Distributions", 1200, 800);
    c2->Divide(3, 2);
    for (int i = 0; i < chargeHistograms.size(); ++i)
    {
        c2->cd(i + 1);
        chargeHistograms[i]->DrawClone();
    }
}