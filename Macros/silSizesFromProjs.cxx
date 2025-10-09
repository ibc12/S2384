#include <map>
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"

void silSizesFromProjs()
{
    std::map<std::string, std::map<int, TH1D *>> hs;

    int nsils{12};
    for (const auto &layer : {"l0", "r0"})
    {
        auto file{std::make_unique<TFile>(TString::Format("../PostAnalysis/Outputs/histos_sp_%s.root", layer))};
        for (int sil = 0; sil < nsils; sil++)
        {
            auto *obj{file->Get(TString::Format("proj%s%d", layer, sil))};
            if (obj)
            {
                auto* h {(TH1D*) obj};
                hs[layer][sil] = h;
                h->SetDirectory(nullptr);
                // Rebin?
                h->Rebin(4);

            }
        }
    }

    // Compute raw heights (no fit)
    auto threshold {4}; // counts threshold
    std::map<std::string, std::map<int, double>> sizes;
    for(auto& [layer, vec] : hs){
        for(auto&[n, h] :  vec){
            auto first {h->FindFirstBinAbove(threshold)};
            auto last {h->FindLastBinAbove(threshold)};
            auto zfirst {h->GetBinCenter(first)};
            auto zlast {h->GetBinCenter(last)};
            auto height {zlast - zfirst};
            sizes[layer][n] = height;
        }
    }
    // Draw
    int canvasIdx{0};
    for (auto &[layer, hsils] : hs)
    {
        // Crear un nuevo canvas para cada histograma
        auto cname = TString::Format("c%d", canvasIdx++);
        auto c = new TCanvas{cname, TString::Format("SP canvas %s", layer.c_str())};
        if (layer == "l0" || layer == "r0")
            c->Divide(3, 4);
        if (layer == "f0")
            c->Divide(4, 3);
        // c0->cd(p);
        int idx{};
        for (auto &[s, h] : hsils)
        {
            c->cd(idx + 1);
            h->SetTitle(TString::Format("%s_%d size: %.2f", layer.c_str(), s, sizes[layer][s]));
            // Draw
            h->Draw("colz");
            idx++;
        }
        canvasIdx++;
    }
}
