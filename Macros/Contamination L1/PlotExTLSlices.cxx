#include "ActCutsManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TColor.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TString.h"
#include "TStyle.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

// Ex in slices (intervals) of an arbitrary column `col` of the dataframe.
// Draws two canvases:
//   cExSlices_<col>  : one pad per slice
//   cExOverlay_<col> : area-normalized overlay of all slices (shape comparison) + Ex vs col (2D)
//
//   nSlices    : number of intervals
//   equalStats : true  -> edges at quantiles of `col` (same number of events per slice)
//                false -> equal-width slices between min and max of `col`
//   userEdges  : if not empty, it overrides the two options above
//   suffix     : appended to the names of canvases / histograms (to be able to call the function several times)
//   note       : text appended to the canvas titles and to the slice histogram titles (e.g. ", RawTL #geq 30 au")
std::vector<double> PlotExSlices(ROOT::RDF::RNode dfin, const std::string& col, const std::string& title,
                                 const std::string& unit, int nSlices, bool equalStats, std::vector<double> userEdges,
                                 const std::string& suffix = "", const std::string& note = "")
{
    const std::string tag {col + suffix};

    // NaN / inf in `col` would break the edges (sort, min/max) and leave every slice empty
    auto df = dfin.Filter([](double x) { return std::isfinite(x); }, {col});

    // ---- Slice edges
    std::vector<double> edges;
    if(!userEdges.empty())
    {
        edges = userEdges;
        std::sort(edges.begin(), edges.end());
        nSlices = static_cast<int>(edges.size()) - 1;
    }
    else
    {
        auto vals {df.Take<double>(col)};
        std::vector<double> v {*vals};
        std::sort(v.begin(), v.end());
        if(v.empty())
        {
            std::cout << "No entries for column " << col << '\n';
            return {};
        }
        if(equalStats)
        {
            for(int i = 0; i <= nSlices; i++)
            {
                std::size_t idx {
                    static_cast<std::size_t>(std::round(static_cast<double>(i) * (v.size() - 1) / nSlices))};
                edges.push_back(v[idx]);
            }
        }
        else
        {
            const double lo {v.front()}, hi {v.back()};
            for(int i = 0; i <= nSlices; i++)
                edges.push_back(lo + (hi - lo) * i / nSlices);
        }
    }
    if(nSlices < 1)
    {
        std::cout << "Need at least 2 edges / 1 slice for " << col << '\n';
        return {};
    }
    // Make the last interval inclusive on the upper edge
    edges.back() = std::nextafter(edges.back(), std::numeric_limits<double>::infinity());
    std::cout << "[" << tag << "] edges:";
    for(const auto& e : edges)
        std::cout << " " << e;
    std::cout << '\n';

    // ---- Histograms (booked lazily; a single pass over the data fills all of them)
    const int nEx {120};
    const double exMin {-5}, exMax {15};

    std::vector<ROOT::RDF::RResultPtr<TH1D>> hs;
    std::vector<TString> labels;
    for(int i = 0; i < nSlices; i++)
    {
        const double lo {edges[i]}, hi {edges[i + 1]};
        auto dfSlice {df.Filter([lo, hi](double x) { return x >= lo && x < hi; }, {col})};
        TString lab {TString::Format("%s #in [%.1f, %.1f) %s", title.c_str(), lo, hi, unit.c_str())};
        labels.push_back(lab);
        hs.push_back(
            dfSlice.Histo1D({TString::Format("hEx_%s_%d", tag.c_str(), i),
                             TString::Format("%s%s;E_{x} [MeV];Counts", lab.Data(), note.c_str()), nEx, exMin, exMax},
                            "Ex"));
    }
    // Reference: Ex vs col (all events)
    auto hExVar {df.Histo2D({TString::Format("hEx_vs_%s", tag.c_str()),
                             TString::Format(";%s [%s];E_{x} [MeV]", title.c_str(), unit.c_str()), 150, edges.front(),
                             edges.back(), nEx, exMin, exMax},
                            col, "Ex")};

    // Entries per slice (first access runs the event loop)
    for(int i = 0; i < nSlices; i++)
        std::cout << "[" << tag << "] slice " << i << ": " << hs[i]->GetEntries()
                  << " entries, integral in Ex range = " << hs[i]->Integral() << '\n';

    const std::vector<int> colors {kBlack,      kRed + 1,  kBlue + 1,   kGreen + 2, kMagenta + 1,
                                   kOrange + 7, kCyan + 2, kViolet + 1, kGray + 2,  kPink + 7};

    // ---- Canvas 1: one pad per slice
    gStyle->SetOptStat(11111);
    const int nCols {static_cast<int>(std::ceil(std::sqrt(static_cast<double>(nSlices))))};
    const int nRows {static_cast<int>(std::ceil(static_cast<double>(nSlices) / nCols))};
    auto* c1 = new TCanvas(TString::Format("cExSlices_%s", tag.c_str()),
                           TString::Format("Ex in %s slices%s", title.c_str(), note.c_str()), 400 * nCols, 350 * nRows);
    c1->Divide(nCols, nRows);
    for(int i = 0; i < nSlices; i++)
    {
        auto* hSl = static_cast<TH1D*>(hs[i]->Clone(TString::Format("hExSlice_%s_%d", tag.c_str(), i)));
        hSl->SetDirectory(nullptr);
        hSl->SetLineWidth(2);
        hSl->SetLineColor(colors[i % colors.size()]);
        c1->cd(i + 1);
        hSl->Draw("hist");
    }
    c1->Modified();
    c1->Update();

    // ---- Canvas 2: overlay (area-normalized) + Ex vs col
    // Los histogramas se copian, se desacoplan de cualquier directorio y se dibujan con Draw()
    // (no DrawClone) para que sigan vivos mientras exista el canvas.
    gStyle->SetOptStat(11111);
    auto* c2 = new TCanvas(TString::Format("cExOverlay_%s", tag.c_str()),
                           TString::Format("Ex overlay and Ex vs %s%s", title.c_str(), note.c_str()), 1600, 700);
    c2->Divide(2, 1);

    c2->cd(1);
    auto* leg = new TLegend(0.55, 0.60, 0.89, 0.89);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    double ymax {0};
    std::vector<TH1D*> norm;
    for(int i = 0; i < nSlices; i++)
    {
        auto* h = static_cast<TH1D*>(hs[i]->Clone(TString::Format("hExNorm_%s_%d", tag.c_str(), i)));
        h->SetDirectory(nullptr);
        if(h->Integral() > 0)
            h->Scale(1. / h->Integral());
        ymax = std::max(ymax, h->GetMaximum());
        h->SetLineColor(colors[i % colors.size()]);
        h->SetLineWidth(2);
        h->SetTitle(";E_{x} [MeV];Normalized counts");
        norm.push_back(h);
        leg->AddEntry(h, TString::Format("%s (N=%d)", labels[i].Data(), static_cast<int>(hs[i]->GetEntries())), "l");
    }
    if(ymax <= 0)
        std::cout << "WARNING: all normalized histograms are empty (ymax = 0) for " << col << '\n';
    for(int i = 0; i < nSlices; i++)
    {
        norm[i]->SetMinimum(0);
        if(std::isfinite(ymax) && ymax > 0)
            norm[i]->SetMaximum(ymax * 1.1);
        norm[i]->Draw(i == 0 ? "hist" : "hist same");
    }
    leg->Draw();

    c2->cd(2);
    auto* h2 = static_cast<TH2D*>(hExVar->Clone(TString::Format("hEx_vs_%s_draw", tag.c_str())));
    h2->SetDirectory(nullptr);
    h2->Draw("colz");

    c1->Modified();
    c1->Update();
    c2->Modified();
    c2->Update();

    return edges;
}

// Ex in slices of the light-particle RawTL and of its azimuthal angle phi (to follow the contamination).
//
//   nSlices / equalStats / userEdges          : slices in RawTL [au]    (equalStats = true  -> same stats per slice)
//   nSlicesPhi / equalStatsPhi / userEdgesPhi : slices in phi light [deg]
//       By default userEdgesPhi is NOT empty, so it overrides nSlicesPhi / equalStatsPhi. The default intervals include
//       [-100, -80) and [80, 100). To use equal-width / equal-stats slices instead, pass userEdgesPhi = {}.
//   Set nSlices = 0 (and no userEdges) to skip the RawTL slices; same for phi (nSlicesPhi = 0 and userEdgesPhi = {}).
//   exContMin / exContMax : Ex window of the contamination [MeV]; the events inside it are overlaid as a scatter
//                           on the PID plots (canvas cPIDscatter)
//   rawTLMin              : minimum RawTL [au] required for the extra Ex-in-phi-slices canvases (suffix
//   _RawTLmin<value>);
//                           pass 0 (or negative) to skip them
//
// The phi intervals are also used for the PID plots (Qtot vs RawTL) of the total data (dftotal).
void PlotExTLSlices(int nSlices = 6, bool equalStats = true, std::vector<double> userEdges = {}, int nSlicesPhi = 6,
                    bool equalStatsPhi = false,
                    std::vector<double> userEdgesPhi = {-180, -135, -100, -80, -45, 0, 45, 80, 100, 135, 180},
                    double exContMin = 5., double exContMax = 6., double rawTLMin = 30.)
{
    std::string beam {"7Li"};
    std::string light {"p"};

    // PID cut (Qtot vs RawTL) for the L1 trigger
    ActRoot::CutsManager<std::string> cuts;
    cuts.ReadCut("l1",
                 TString::Format("../../PostAnalysis/Cuts/pid_%s_l1_%s.root", light.c_str(), beam.c_str()).Data());

    ROOT::RDataFrame dftotal_ini {
        "PreProcessed_Tree", TString::Format("../../PostAnalysis/Outputs/tree_preprocess_F_%s.root", beam.c_str())};
    auto dftotal = dftotal_ini.Filter([](ActRoot::ModularData& m) { return m.Get("GATCONF") == 8; }, {"ModularData"});
    std::cout << dftotal.Count().GetValue() << " events in the tree" << '\n';

    ROOT::RDataFrame dfini {"Final_Tree", TString::Format("../../PostAnalysis/Outputs/tree_ex_F_%s_d_%s_filtered.root",
                                                          beam.c_str(), light.c_str())};

    auto df = dfini
                  .Define("RawTL", [](const ActRoot::MergerData& d) { return static_cast<double>(d.fLight.fRawTL); },
                          {"MergerData"})
                  // NOTE: assumed member name (in degrees, like fThetaLight). Adjust if it is called differently
                  .Define("PhiL", [](const ActRoot::MergerData& d) { return static_cast<double>(d.fPhiLight); },
                          {"MergerData"})
                  .Define("Qtot", [](const ActRoot::MergerData& d) { return static_cast<double>(d.fLight.fQtotal); },
                          {"MergerData"})
                  .Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"});

    // ---- Ex in slices of RawTL
    if(nSlices > 0 || !userEdges.empty())
        PlotExSlices(df, "RawTL", "RawTL", "au", nSlices, equalStats, userEdges);

    // ---- Ex in slices of phi of the light particle
    std::vector<double> phiEdges;
    if(nSlicesPhi > 0 || !userEdgesPhi.empty())
        phiEdges = PlotExSlices(df, "PhiL", "#phi_{L}", "deg", nSlicesPhi, equalStatsPhi, userEdgesPhi);

    // ---- Same phi slices, but only with events with RawTL >= rawTLMin
    if(rawTLMin > 0 && (nSlicesPhi > 0 || !userEdgesPhi.empty()))
    {
        auto dfMin {df.Filter([rawTLMin](double raw) { return raw >= rawTLMin; }, {"RawTL"})};
        PlotExSlices(dfMin, "PhiL", "#phi_{L}", "deg", nSlicesPhi, equalStatsPhi, userEdgesPhi,
                     TString::Format("_RawTLmin%g", rawTLMin).Data(),
                     TString::Format(", RawTL #geq %g au", rawTLMin).Data());
    }

    // ---- Canvas 3: Qtotal vs TL and Qtotal vs RawTL (dftotal, preprocessed tree, GATCONF == 8)
    auto dfq = dftotal
                   .Define("Qtot", [](const ActRoot::MergerData& d) { return static_cast<double>(d.fLight.fQtotal); },
                           {"MergerData"})
                   .Define("TLq", [](const ActRoot::MergerData& d) { return static_cast<double>(d.fLight.fTL); },
                           {"MergerData"})
                   .Define("RawTL", [](const ActRoot::MergerData& d) { return static_cast<double>(d.fLight.fRawTL); },
                           {"MergerData"})
                   .Define("PhiL", [](const ActRoot::MergerData& d) { return static_cast<double>(d.fPhiLight); },
                           {"MergerData"});

    // Ranges from the data (extra pass)
    auto qMin {dfq.Min<double>("Qtot")};
    auto qMax {dfq.Max<double>("Qtot")};
    auto tlMin {dfq.Min<double>("TLq")};
    auto tlMax {dfq.Max<double>("TLq")};
    auto rawMin {dfq.Min<double>("RawTL")};
    auto rawMax {dfq.Max<double>("RawTL")};
    const double q0 {*qMin}, q1 {*qMax * 1.05};
    const double t0 {*tlMin}, t1 {*tlMax * 1.05};
    const double r0 {*rawMin}, r1 {*rawMax * 1.05};

    const int nQ {150}, nT {150};
    auto hQTL {dfq.Histo2D({"hQTLtot", ";TL [mm];Q_{tot} [arb. u.]", nT, t0, t1, 2000, 0, 3e5}, "TLq", "Qtot")};
    auto hQRaw {dfq.Histo2D({"hQRawTL", ";RawTL [au];Q_{tot} [arb. u.]", nT, r0, r1, 2000, 0, 3e5}, "RawTL", "Qtot")};

    gStyle->SetOptStat(11111);
    auto* c3 = new TCanvas("cQTL", "Qtotal vs TL and RawTL", 1600, 700);
    c3->Divide(2, 1);
    c3->cd(1);
    auto* h3a = static_cast<TH2D*>(hQTL->Clone("hQTL_draw"));
    h3a->SetDirectory(nullptr);
    h3a->Draw("colz");
    c3->cd(2);
    auto* h3b = static_cast<TH2D*>(hQRaw->Clone("hQRawTL_draw"));
    h3b->SetDirectory(nullptr);
    h3b->Draw("colz");
    cuts.DrawCut("l1");

    // ---- Canvas 4: PID (Qtot vs RawTL) for each phi interval, total data (dftotal)
    if(phiEdges.size() >= 2)
    {
        const int nPhi {static_cast<int>(phiEdges.size()) - 1};
        const int nQslice {400}; // Q bins per histogram (range 0 - 3e5, as in canvas 3)
        std::vector<ROOT::RDF::RResultPtr<TH2D>> hPID;
        for(int i = 0; i < nPhi; i++)
        {
            const double lo {phiEdges[i]}, hi {phiEdges[i + 1]};
            auto dfPhi {dfq.Filter([lo, hi](double phi) { return phi >= lo && phi < hi; }, {"PhiL"})};
            TString lab {TString::Format("#phi_{L} #in [%g, %g) deg", lo, hi)};
            hPID.push_back(dfPhi.Histo2D({TString::Format("hPID_phi_%d", i),
                                          TString::Format("%s;RawTL [au];Q_{tot} [arb. u.]", lab.Data()), nT, r0, r1,
                                          nQslice, 0, 3e5},
                                         "RawTL", "Qtot"));
        }
        // First access runs the event loop for all the histograms
        for(int i = 0; i < nPhi; i++)
            std::cout << "[PID] phi slice " << i << ": " << hPID[i]->GetEntries() << " entries" << '\n';

        gStyle->SetOptStat(11111);
        const int nColsP {static_cast<int>(std::ceil(std::sqrt(static_cast<double>(nPhi))))};
        const int nRowsP {static_cast<int>(std::ceil(static_cast<double>(nPhi) / nColsP))};
        auto* c4 = new TCanvas("cPIDphi", "Qtot vs RawTL in phi intervals (dftotal)", 450 * nColsP, 380 * nRowsP);
        c4->Divide(nColsP, nRowsP);
        for(int i = 0; i < nPhi; i++)
        {
            auto* hp = static_cast<TH2D*>(hPID[i]->Clone(TString::Format("hPID_phi_%d_draw", i)));
            hp->SetDirectory(nullptr);
            c4->cd(i + 1);
            hp->Draw("colz");
            cuts.DrawCut("l1");
            // c4->cd(i + 1)->SetLogz();
        }
        c4->Modified();
        c4->Update();

        // ---- Canvas 5: same PID as canvas 4 (+ all phi in the first pad) with a scatter on top of the events with
        // exContMin <= Ex < exContMax (contamination region), to see where they sit with respect to the PID cut.
        // Ex only exists in the final tree (df, L1 events), so the points come from df and the PID from dftotal.
        auto dfCont {
            df.Filter([exContMin, exContMax](double ex) { return ex >= exContMin && ex < exContMax; }, {"Ex"})};
        auto contRaw {dfCont.Take<double>("RawTL")};
        auto contQ {dfCont.Take<double>("Qtot")};
        auto contPhi {dfCont.Take<double>("PhiL")};
        // First access runs the event loop for the three columns
        const std::vector<double> vRaw {*contRaw};
        const std::vector<double> vQ {*contQ};
        const std::vector<double> vPhi {*contPhi};
        std::cout << "[Scatter] " << vRaw.size() << " events with Ex in [" << exContMin << ", " << exContMax << ") MeV"
                  << '\n';

        // Graph with the contamination events inside a phi interval (or all of them)
        auto makeGraph = [&](bool useAll, double phiLo, double phiHi)
        {
            auto* g = new TGraph();
            for(std::size_t k = 0; k < vRaw.size(); k++)
            {
                if(!useAll && !(vPhi[k] >= phiLo && vPhi[k] < phiHi))
                    continue;
                g->SetPoint(g->GetN(), vRaw[k], vQ[k]);
            }
            g->SetMarkerStyle(20);
            g->SetMarkerSize(0.5);
            g->SetMarkerColor(kRed + 1);
            return g;
        };
        auto drawScatter = [&](TGraph* g)
        {
            if(g->GetN() == 0)
                return;
            g->Draw("P same");
            auto* lg = new TLegend(0.12, 0.80, 0.70, 0.88);
            lg->SetBorderSize(0);
            lg->SetFillStyle(0);
            lg->AddEntry(g, TString::Format("E_{x} #in [%g, %g) MeV (N=%d)", exContMin, exContMax, g->GetN()), "p");
            lg->Draw();
        };

        const int nPads {nPhi + 1};
        const int nColsS {static_cast<int>(std::ceil(std::sqrt(static_cast<double>(nPads))))};
        const int nRowsS {static_cast<int>(std::ceil(static_cast<double>(nPads) / nColsS))};
        auto* c5 =
            new TCanvas("cPIDscatter", "PID + scatter of the contamination (Ex window)", 450 * nColsS, 380 * nRowsS);
        c5->Divide(nColsS, nRowsS);

        // Pad 1: all phi
        c5->cd(1);
        auto* hAll = static_cast<TH2D*>(hQRaw->Clone("hQRawTL_scatter"));
        hAll->SetDirectory(nullptr);
        hAll->SetTitle("All #phi_{L};RawTL [au];Q_{tot} [arb. u.]");
        hAll->Draw("colz");
        cuts.DrawCut("l1");
        drawScatter(makeGraph(true, 0., 0.));

        // Pads 2...: one per phi interval
        for(int i = 0; i < nPhi; i++)
        {
            c5->cd(i + 2);
            auto* hp = static_cast<TH2D*>(hPID[i]->Clone(TString::Format("hPID_phi_%d_scatter", i)));
            hp->SetDirectory(nullptr);
            hp->Draw("colz");
            cuts.DrawCut("l1");
            drawScatter(makeGraph(false, phiEdges[i], phiEdges[i + 1]));
        }
        c5->Modified();
        c5->Update();
    }
}