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
    std::vector<double> maxCharges {};
    for (int run = 124; run <= 128; ++run)
    {
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
        file.GetObject("GETTree", tree); // cambia "tree" por el nombre real
        if (!tree) {
            std::cerr << "No se encontró TTree en " << filename << "\n";
            continue;
        }

        // Conectar branch de TPCData
        ActRoot::TPCData *tpcData = nullptr;
        tree->SetBranchAddress("TPCData", &tpcData);

        double sumCharge = 0;
        Long64_t nEvents = tree->GetEntries();

        // Recorrer eventos
        for (Long64_t i = 0; i < nEvents; ++i)
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
            sumCharge += maxCharge;
        }

        double meanCharge = (nEvents > 0) ? sumCharge / nEvents : 0;
        maxCharges.push_back(meanCharge);
        std::cout << "Run " << run << " -> media maxCharge = " << meanCharge << "\n";
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
    graph->Fit("pol2", "0Q"); // Fit with a polynomial of degree 2
    TF1 *fitFunc = graph->GetFunction("pol2");
    // Plot graph
    TCanvas *c1 = new TCanvas("c1", "Gain vs Mesh Voltage", 800, 600);
    graph->DrawClone("AP");
    fitFunc->SetLineColor(kRed);
    fitFunc->Draw("same");
}