#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TSpectrum.h"
#include "TGraph.h"
#include "TF1.h"
#include "TStyle.h"
#include <fstream>
#include <iostream>

std::vector<double> GetPeaks(TH2D *h, int c)
{
    auto *proj{h->ProjectionY("proj", c, c)};

    auto *spe{new TSpectrum(11)}; // max nums of peaks, is overstimated
    int nPeaks{spe->Search(proj, 2, "nodraw", 0.05)};

    double *xPositions{spe->GetPositionX()}; // We have to sort them

    std::sort(xPositions, xPositions + nPeaks); // std::sort needs a pointer for the begin and end of the vector

    std::vector<double> peaks(xPositions, xPositions + nPeaks); // convert the pointer to a vector to easier use

    double width{50.};
    std::vector<double> meanPeak;

    for (auto &peak : peaks)
    {
        double xmin{peak - width};
        double xmax{peak + width};

        proj->Fit("gaus", "0Q", "", xmin, xmax);
        meanPeak.push_back(proj->GetFunction("gaus")->GetParameter("Mean"));
    }

    return meanPeak;
}

void FillGraph(TGraph *graph, const std::vector<double> &x, const std::vector<double> &y, std::ofstream &streamer)
{
    if (x.size() != y.size())
    {
        streamer << 0 << " " << 0 << " " << 0 << std::endl;
        return;
    }
    for (int p = 0; p < x.size(); p++)
    {
        graph->SetPoint(p, x[p], y[p]);
    }

    graph->Fit("pol2", "0Q");
    TF1 *f{graph->GetFunction("pol2")};
    if (!f)
    {
        std::cout << "Func doesn't exist" << '\n';
        streamer << 0 << " " << 0 << " " << 0 << std::endl;
    }
    else
    {
        streamer << f->GetParameter(0) << " " << f->GetParameter(1) << " " << f->GetParameter(2) << std::endl;
    }
}

void DoGainMatching()
{
    auto *file{new TFile{"./Inputs/gain.root"}};
    auto *h{file->Get<TH2D>("h")};

    int channels{17408};

    std::vector<std::vector<double>> peaks;

    for (int c = 0; c < channels; c++)
    {
        peaks.push_back(GetPeaks(h, c));
    }

    int channelRef{4196};

    // Now we have to do the fit Q vs Q ref, we use a TGraph, for filling it need a for loop

    std::ofstream streamer{"./Outputs/gain_matching_v0.dat"};

    std::vector<TGraph *> gs;
    for (const auto &peak : peaks)
    {
        TGraph *graph{new TGraph()};
        FillGraph(graph, peak, peaks[channelRef-1], streamer); //channelRef-1 takes into acount that the bins start in 1 but channels in 0
        graph->SetMarkerStyle(24);
        gs.push_back(graph);
    }
    streamer.close();

    gStyle->SetOptFit(true);
    auto *c{new TCanvas("c")};
    int chosen{4196};
    gs[chosen]->Draw("ap");
    for (auto *ptr : *(gs[chosen]->GetListOfFunctions()))
        if (ptr)
            ptr->Draw("same");
}