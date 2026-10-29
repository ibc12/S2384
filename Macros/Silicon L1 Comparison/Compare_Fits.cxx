// Compare_Fits.C
// Compara los resultados de dos fits guardados con Fitters::RunFit (Runner::Write + SaveGlobalFit).
// Solo se comparan los estados que existen en AMBOS ficheros (mismo nombre: g0, g1, v0, ...).
//
// Uso:
//   root -l 'Compare_Fits.C("./Outputs/fit_A.root", "./Outputs/fit_B.root", "Fit A", "Fit B")'
//
// Se espera en cada fichero:
//   ParNames    -> std::vector<std::string>  (nombre del parametro; la clave del estado es lo anterior al primer '_')
//   FitResult   -> TFitResult
//   HistoEx     -> TH1D   (espectro de Ex)
//   GraphGlobal -> TGraph (ajuste global)

#include "TCanvas.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TLegend.h"
#include "TString.h"
#include "TStyle.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace CompareFits
{
struct FitData
{
    std::string fName;
    std::string fLabel;
    std::map<std::string, std::vector<double>> fVals; // key -> parameters of that state
    std::map<std::string, std::vector<double>> fUncs;
    double fMinFcn {-1};
    int fNdf {-1};
    TH1D* fHEx {};
    TGraph* fGlobal {};
    bool fOk {false};
};

// ---------- Helpers ----------
int Group(const std::string& k)
{
    if(k.rfind("ps", 0) == 0)
        return 2;
    if(!k.empty() && k[0] == 'g')
        return 0;
    if(!k.empty() && k[0] == 'v')
        return 1;
    return 3;
}

int Number(const std::string& k)
{
    auto p {k.find_first_of("0123456789")};
    return p == std::string::npos ? -1 : std::atoi(k.c_str() + p);
}

bool KeyLess(const std::string& a, const std::string& b)
{
    if(Group(a) != Group(b))
        return Group(a) < Group(b);
    if(Number(a) != Number(b))
        return Number(a) < Number(b);
    return a < b;
}

// Order of parameters: gauss {N, mean, sigma}; voigt {N, mean, sigma, gamma}; ps {scale}
std::string ParName(const std::string& key, int i)
{
    if(Group(key) == 2)
        return "Scale";
    static const std::vector<std::string> names {"Amplitude", "Mean", "Sigma", "Gamma"};
    return i < (int)names.size() ? names[i] : ("p" + std::to_string(i));
}

// ---------- Read ----------
FitData Load(const std::string& file, const std::string& label)
{
    FitData d;
    d.fName = file;
    d.fLabel = label;
    std::unique_ptr<TFile> f {TFile::Open(file.c_str())};
    if(!f || f->IsZombie())
    {
        std::cout << "Compare_Fits: cannot open " << file << '\n';
        return d;
    }
    auto* names {f->Get<std::vector<std::string>>("ParNames")};
    auto* res {f->Get<TFitResult>("FitResult")};
    if(!names || !res)
    {
        std::cout << "Compare_Fits: ParNames or FitResult missing in " << file << '\n';
        return d;
    }
    int n {std::min((int)names->size(), (int)res->NPar())};
    for(int i = 0; i < n; i++)
    {
        auto it {names->at(i).find_first_of("_")};
        auto key {names->at(i).substr(0, it)};
        d.fVals[key].push_back(res->Parameter(i));
        d.fUncs[key].push_back(res->Error(i));
    }
    d.fMinFcn = res->MinFcnValue();
    d.fNdf = res->Ndf();
    // Histogram and global fit (survive file closing)
    d.fHEx = f->Get<TH1D>("HistoEx");
    if(d.fHEx)
        d.fHEx->SetDirectory(nullptr);
    d.fGlobal = f->Get<TGraph>("GraphGlobal");
    d.fOk = true;
    return d;
}

std::vector<std::string> CommonKeys(const FitData& a, const FitData& b)
{
    std::vector<std::string> keys;
    for(const auto& [k, v] : a.fVals)
    {
        if(!b.fVals.count(k))
            continue;
        if(v.size() != b.fVals.at(k).size())
        {
            std::cout << "Compare_Fits: '" << k << "' has different number of parameters in both fits -> skipped\n";
            continue;
        }
        keys.push_back(k);
    }
    std::sort(keys.begin(), keys.end(), KeyLess);
    return keys;
}

void PrintOnlyIn(const FitData& a, const FitData& b)
{
    std::vector<std::string> only;
    for(const auto& [k, v] : a.fVals)
        if(!b.fVals.count(k))
            only.push_back(k);
    std::sort(only.begin(), only.end(), KeyLess);
    if(only.empty())
        return;
    std::cout << "States only in '" << a.fLabel << "' (not compared): ";
    for(const auto& k : only)
        std::cout << k << ' ';
    std::cout << '\n';
}

// ---------- Table ----------
void PrintTable(const FitData& a, const FitData& b, const std::vector<std::string>& keys)
{
    std::printf("\n%s: %s  (FCN/ndf = %.2f / %d)\n", a.fLabel.c_str(), a.fName.c_str(), a.fMinFcn, a.fNdf);
    std::printf("%s: %s  (FCN/ndf = %.2f / %d)\n", b.fLabel.c_str(), b.fName.c_str(), b.fMinFcn, b.fNdf);
    std::printf("\n%-6s %-10s %-20s %-20s %-12s %-8s\n", "State", "Param", a.fLabel.c_str(), b.fLabel.c_str(),
                "Diff(B-A)", "nSigma");
    std::printf("-------------------------------------------------------------------------------------\n");
    for(const auto& k : keys)
    {
        const auto& va {a.fVals.at(k)};
        const auto& vb {b.fVals.at(k)};
        const auto& ua {a.fUncs.at(k)};
        const auto& ub {b.fUncs.at(k)};
        for(size_t i = 0; i < va.size(); i++)
        {
            double diff {vb[i] - va[i]};
            double den {std::sqrt(ua[i] * ua[i] + ub[i] * ub[i])};
            TString sa {TString::Format("%.4g +- %.2g", va[i], ua[i])};
            TString sb {TString::Format("%.4g +- %.2g", vb[i], ub[i])};
            TString ns {den > 0 ? TString::Format("%.2f", diff / den) : TString("--")};
            const char* flag {(den > 0 && std::abs(diff / den) > 2) ? " *" : ""};
            std::printf("%-6s %-10s %-20s %-20s %-12.4g %-8s%s\n", (i == 0 ? k.c_str() : ""), ParName(k, i).c_str(),
                        sa.Data(), sb.Data(), diff, ns.Data(), flag);
        }
    }
    std::printf("\n(* = diferencia > 2 sigma; '--' = ambos parametros fijos / sin error)\n\n");
}

// ---------- Plots ----------
void DrawGlobal(const FitData& a, const FitData& b)
{
    if(!a.fHEx || !b.fHEx || !a.fGlobal || !b.fGlobal)
    {
        std::cout << "Compare_Fits: HistoEx / GraphGlobal missing, global overlay skipped\n";
        return;
    }
    auto* c {new TCanvas {"cCmpGlobal", "Global fit comparison", 1000, 700}};
    auto* h1 {(TH1D*)a.fHEx->Clone("hCmpA")};
    auto* h2 {(TH1D*)b.fHEx->Clone("hCmpB")};
    h1->SetDirectory(nullptr);
    h2->SetDirectory(nullptr);
    for(auto* h : {h1, h2})
    {
        h->SetStats(false);
        h->SetTitle("");
        h->SetLineWidth(2);
    }
    h1->SetLineColor(kBlack);
    h2->SetLineColor(kGray + 1);
    // Draw first the histogram with the wider axis so that nothing is cut
    auto width {[](TH1D* h) { return h->GetXaxis()->GetXmax() - h->GetXaxis()->GetXmin(); }};
    auto* first {width(h1) >= width(h2) ? h1 : h2};
    auto* second {first == h1 ? h2 : h1};
    double ymax {1.2 * std::max(h1->GetMaximum(), h2->GetMaximum())};
    first->SetMaximum(ymax);
    first->Draw("e");
    second->Draw("e same");

    auto* g1 {(TGraph*)a.fGlobal->Clone()};
    auto* g2 {(TGraph*)b.fGlobal->Clone()};
    g1->SetLineColor(kAzure + 2);
    g2->SetLineColor(kRed + 1);
    g1->SetLineWidth(3);
    g2->SetLineWidth(3);
    g2->SetLineStyle(7);
    g1->Draw("l same");
    g2->Draw("l same");

    auto* leg {new TLegend {0.5, 0.65, 0.9, 0.9}};
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h1, TString::Format("Data %s", a.fLabel.c_str()), "le");
    leg->AddEntry(h2, TString::Format("Data %s", b.fLabel.c_str()), "le");
    leg->AddEntry(g1, TString::Format("%s (FCN/ndf = %.1f/%d)", a.fLabel.c_str(), a.fMinFcn, a.fNdf), "l");
    leg->AddEntry(g2, TString::Format("%s (FCN/ndf = %.1f/%d)", b.fLabel.c_str(), b.fMinFcn, b.fNdf), "l");
    leg->Draw();
    c->Update();
}

void DrawParameters(const FitData& a, const FitData& b, const std::vector<std::string>& keys)
{
    // Only states with peaks (no phase-space scales)
    std::vector<std::string> states;
    for(const auto& k : keys)
        if(Group(k) != 2)
            states.push_back(k);
    if(states.empty())
        return;

    auto* c {new TCanvas {"cCmpPars", "Parameter comparison", 1400, 900}};
    c->Divide(2, 2);
    const std::vector<std::string> parTitles {"Amplitude", "Mean (MeV)", "Sigma (MeV)", "Gamma (MeV)"};
    const int colA {kAzure + 2};
    const int colB {kRed + 1};
    const double shift {0.1};

    for(int ip = 0; ip < 4; ip++)
    {
        c->cd(ip + 1);
        auto* ga {new TGraphErrors};
        auto* gb {new TGraphErrors};
        double ymin {1e300};
        double ymax {-1e300};
        for(size_t is = 0; is < states.size(); is++)
        {
            const auto& k {states[is]};
            const auto& va {a.fVals.at(k)};
            const auto& vb {b.fVals.at(k)};
            if((int)va.size() <= ip)
                continue; // e.g. no gamma for gaussians
            const auto& ua {a.fUncs.at(k)};
            const auto& ub {b.fUncs.at(k)};
            double x {is + 0.5};
            int n {ga->GetN()};
            ga->SetPoint(n, x - shift, va[ip]);
            ga->SetPointError(n, 0, ua[ip]);
            gb->SetPoint(n, x + shift, vb[ip]);
            gb->SetPointError(n, 0, ub[ip]);
            ymin = std::min({ymin, va[ip] - ua[ip], vb[ip] - ub[ip]});
            ymax = std::max({ymax, va[ip] + ua[ip], vb[ip] + ub[ip]});
        }
        if(ga->GetN() == 0)
            continue;
        double pad {0.15 * (ymax - ymin)};
        if(pad <= 0)
            pad = 0.1 * std::max(std::abs(ymax), 1e-3);
        // Frame with state labels
        auto* frame {new TH1F {TString::Format("frameCmp%d", ip), TString::Format(";;%s", parTitles[ip].c_str()),
                               (int)states.size(), 0, (double)states.size()}};
        frame->SetDirectory(nullptr);
        frame->SetStats(false);
        for(size_t is = 0; is < states.size(); is++)
            frame->GetXaxis()->SetBinLabel(is + 1, states[is].c_str());
        frame->GetXaxis()->SetLabelSize(0.06);
        frame->SetMinimum(ymin - pad);
        frame->SetMaximum(ymax + pad);
        frame->Draw("axis");
        for(auto* g : {ga, gb})
        {
            g->SetMarkerStyle(20);
            g->SetMarkerSize(1.1);
            g->SetLineWidth(2);
        }
        ga->SetMarkerColor(colA);
        ga->SetLineColor(colA);
        gb->SetMarkerColor(colB);
        gb->SetLineColor(colB);
        ga->Draw("pz same");
        gb->Draw("pz same");
        if(ip == 0)
        {
            auto* leg {new TLegend {0.65, 0.78, 0.92, 0.92}};
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(ga, a.fLabel.c_str(), "lp");
            leg->AddEntry(gb, b.fLabel.c_str(), "lp");
            leg->Draw();
        }
    }
    c->Update();
}
} // namespace CompareFits

// This is a macro to compare the result of the Ex spectrum fitting between the silicon and L1 spectra

void Compare_Fits(const std::string& file1 = "../../Fits/7Li_dp/Outputs/fit_lDependent.root",
                  const std::string& file2 = "../../Fits/7Li_dp_testL1/Outputs/fit_lDependent.root", const std::string& label1 = "Fit silicon",
                  const std::string& label2 = "Fit L1")
{
    using namespace CompareFits;
    gStyle->SetOptStat(0);

    auto a {Load(file1, label1)};
    auto b {Load(file2, label2)};
    if(!a.fOk || !b.fOk)
        return;

    auto keys {CommonKeys(a, b)};
    if(keys.empty())
    {
        std::cout << "Compare_Fits: no common states between both fits\n";
        return;
    }
    PrintOnlyIn(a, b);
    PrintOnlyIn(b, a);
    PrintTable(a, b, keys);
    DrawGlobal(a, b);
    DrawParameters(a, b, keys);
}