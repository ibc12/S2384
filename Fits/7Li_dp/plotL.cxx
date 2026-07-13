#include "TCanvas.h"
#include "TGraph.h"

#include <fstream>
#include <iostream>
#include <vector>

void plotL()
{

    // Read.dat from freasco output, first and thrird column

    std::ifstream infile("./Inputs/gs_ADWA/fort.56");
    std::vector<double> l, sigma;
    double l_val, sigma_val, sigma_1, sigma_2;
    while(infile >> l_val >> sigma_1 >> sigma_2 >> sigma_val)
    {
        // skip first row
        l.push_back(l_val);
        sigma.push_back(sigma_val);
    }
    std::cout << "Read " << l.size() << " points from fort.56\n";
    for(const auto& val : l)
    {
        std::cout << val << " ";
    }

    // Create a TGraph
    TGraph* graph = new TGraph(l.size(), &l[0], &sigma[0]);
    graph->SetTitle("Cross Section vs L;L;Cross Section");
    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(kBlue);

    // Draw the graph
    TCanvas* c1 = new TCanvas("c1", "Cross Section vs L", 800, 600);
    graph->Draw("AP");
}