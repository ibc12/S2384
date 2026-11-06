// Toy: construye un histograma de Ex con los estados definidos en el Interface
// (con sus parametros INICIALES: amplitud, media, sigma, gamma) mas el espacio fasico
// (hPS escalado por la amplitud inicial de ps0), y lo ajusta con RunFit.
//
//   root -l 'Toy_lDependent.C'                  -> toy con fluctuacion Poisson, fit desde los valores iniciales
//   root -l 'Toy_lDependent.C(false)'           -> histograma exacto (sin fluctuar)
//   root -l 'Toy_lDependent.C(true, 1, 0.2)'    -> seed 1, y el fit arranca con amplitudes/medias perturbadas un ~20%

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TString.h"

#include "FitInterface.h"
#include "FitModel.h"
#include "FitRunner.h"
#include "FitUtils.h"
#include "Interpolators.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "../../PrettyStyle.C"
#include "../Histos.h"

// Valor esperado por bin: promedio del modelo en nsub subdivisiones del bin
// (igual que SetUseDivisions(true) en el objetivo del fit)
TH1D* ToyExpected(Fitters::Model& model, std::vector<double> pars, const TH1D& tmpl, const char* name, double exmin,
                  double exmax, int nsub = 10)
{
    auto* h {static_cast<TH1D*>(tmpl.Clone(name))};
    h->SetDirectory(nullptr);
    h->Reset();
    h->Sumw2(false);
    // Los Voigt con AddBWL se evaluan con TF1Convolution (BW(l) x Gauss), que solo existe tras
    // TriggerConvolution(). Hay que llamarlo con ESTOS parametros antes de evaluar el modelo.
    model.TriggerConvolution(pars.data(), exmin, exmax);
    std::cout << "[Toy] TriggerConvolution llamado para " << name << '\n';
    for(int b = 1; b <= h->GetNbinsX(); b++)
    {
        if(h->GetBinCenter(b) < exmin || h->GetBinCenter(b) > exmax)
            continue;
        double lo {h->GetBinLowEdge(b)};
        double w {h->GetBinWidth(b)};
        double mu {};
        for(int k = 0; k < nsub; k++)
        {
            double x {lo + (k + 0.5) * w / nsub};
            mu += model(&x, pars.data());
        }
        h->SetBinContent(b, mu / nsub);
    }
    return h;
}

void Toy_lDependent(bool fluctuate = true, int seed = 12345, double jitter = 0.)
{
    PrettyStyle(false, false);
    TRandom3 rng {static_cast<unsigned int>(seed)};

    // Fitting range
    double exmin {-1};
    double exmax {10};

    // ---------------- Phase space (igual que Fit_lDependent) ----------------
    ROOT::RDataFrame phase {"SimulationTTree", "../../Simulation/Outputs/7Li/2H_1H_TRIUMF_Eex_0.000_nPS_1_pPS_0.root"};
    auto hPS {phase.Histo1D(S2384Fit::Exdp_7Li, "Eex_side", "weight")};
    hPS->SetNameTitle("hPS", "1n PS");

    // ---------------- Sigmas e Interface (copia de Fit_lDependent) ----------------
    Interpolators::Sigmas sigmas;
    sigmas.Read("../../Simulation/Outputs/7Li/test_ang_straggling/sigmas_7Li_2H_1H_1-4AngStr.root");

    Fitters::Interface inter;
    double sigma {0.135};
    inter.AddState("g0", {100, 0, sigma});
    inter.AddState("g1", {30, 0.98, sigma});
    inter.AddState("v0", {60, 2.2, sigma, 0.1});
    inter.AddState("v1", {5, 2.85, sigma, 0.4});
    inter.AddState("v2", {20, 3.3, sigma, 1});
    inter.AddState("v3", {15, 3.9, sigma, 0.65});
    inter.AddState("v4", {15, 5.2, sigma, 0.65});
    inter.AddState("v5", {15, 5.7, sigma, 1});
    inter.AddState("v6", {10, 6.3, sigma, 0.035});
    inter.AddState("v7", {5, 7., sigma, 0.4});
    inter.AddState("v8", {5, 7.6, sigma, 0.4});
    inter.AddState("ps0", {0.0005});
    inter.EndAddingStates();
    inter.EvalSigma(sigmas.GetGraph());
    inter.SetFixAll(2, true); // fix all sigmas
    inter.SetBounds("v0", 3, {0.005, 0.3});
    inter.SetBounds("v1", 3, {0.1, 0.8});
    inter.SetBounds("v2", 3, {0.5, 1.5});
    inter.SetBounds("v3", 3, {0.01, 1});
    inter.SetBounds("v4", 3, {0.3, 1});
    inter.SetBounds("v5", 3, {0.5, 1.5});
    inter.SetBounds("v6", 3, {0.01, 0.1});
    inter.SetBounds("v7", 3, {0.1, 0.8});
    inter.SetBounds("v8", 3, {0.1, 0.8});
    inter.SetOffsetMeanBounds(0.25);

    // ---------------- Model (igual que Fit_lDependent) ----------------
    Fitters::Model model {inter.GetNGauss(), inter.GetNVoigt(), {*hPS}};
    double R {(std::pow(7, 1. / 3.) + std::pow(1, 1. / 3.)) * 1.25};
    double mu {7 * 1. / (7 + 1) * 931.5};
    model.AddBWL(0, 1, 2.03262, mu, R);
    model.AddBWL(1, 1, 2.03262, mu, R);
    model.AddBWL(2, 1, 2.03262, mu, R);
    model.AddBWL(3, 1, 2.03262, mu, R);
    model.AddBWL(4, 2, 2.03262, mu, R);
    model.AddBWL(5, 2, 2.03262, mu, R);
    model.AddBWL(6, 2, 2.03262, mu, R);
    model.AddBWL(7, 1, 2.03262, mu, R);
    model.AddBWL(8, 1, 2.03262, mu, R);

    // ---------------- Vector plano de parametros iniciales ----------------
    // Se asume que el orden del modelo es el de los estados (g0, g1, v0..v8, ps0),
    // con los parametros de cada uno seguidos. Se comprueba contra ParNames tras el fit.
    auto init {inter.GetInitial()};
    const auto keys {inter.GetKeys()};
    std::vector<double> pars;
    std::vector<std::string> keyOfPar;    // key a la que pertenece cada parametro
    std::map<std::string, size_t> ampIdx; // posicion de la amplitud de cada key
    for(const auto& key : keys)
    {
        const auto& v {init.at(key)};
        ampIdx[key] = pars.size();
        for(auto x : v)
        {
            pars.push_back(x);
            keyOfPar.push_back(key);
        }
    }
    std::cout << "Parametros iniciales: " << pars.size() << " | model.NPar(): " << model.NPar() << '\n';
    if(pars.size() != static_cast<size_t>(model.NPar()))
    {
        std::cout << "!! El numero de parametros no coincide: revisa el orden/numero de parametros del modelo" << '\n';
        return;
    }

    // ---------------- Histograma esperado + toy ----------------
    // Plantilla de binning: la del phase space (mismo Exdp_7Li que los datos)
    auto* hExp {ToyExpected(model, pars, *hPS, "hExpected", exmin, exmax)};

    auto* hToy {static_cast<TH1D*>(hExp->Clone("hToy"))};
    hToy->SetDirectory(nullptr);
    hToy->SetTitle("Toy;E_{x} [MeV];Counts");
    if(fluctuate)
    {
        for(int b = 1; b <= hToy->GetNbinsX(); b++)
            hToy->SetBinContent(b, rng.Poisson(hExp->GetBinContent(b)));
    }
    std::cout << "Integral esperado = " << hExp->Integral() << " | toy = " << hToy->Integral() << '\n';

    // Componentes individuales (solo para dibujar y comprobar): amplitud de las demas a 0
    auto* c0 {new TCanvas {"cToyGen", "Toy generated"}};
    hToy->SetStats(false);
    hToy->DrawClone("e");
    auto* leg {new TLegend {0.55, 0.55, 0.9, 0.9}};
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetNColumns(3);
    int icol {};
    for(const auto& key : keys)
    {
        auto p {pars};
        for(const auto& [other, idx] : ampIdx)
            if(other != key)
                p[idx] = 0;
        auto* hc {ToyExpected(model, p, *hPS, ("hComp_" + key).c_str(), exmin, exmax)};
        hc->SetLineColor(2 + icol % 8);
        hc->SetLineWidth(2);
        hc->DrawClone("hist same");
        leg->AddEntry(hc->Clone(), key.c_str(), "l");
        icol++;
    }
    leg->Draw();
    std::cout << "Counts del phase space (ps0): "
              << ToyExpected(
                     model,
                     [&]()
                     {
                         auto p {pars};
                         for(const auto& [other, idx] : ampIdx)
                             if(other != "ps0")
                                 p[idx] = 0;
                         return p;
                     }(),
                     *hPS, "hPSonly", exmin, exmax)
                     ->Integral()
              << '\n';

    // ---------------- Fit con RunFit ----------------
    auto fitInit {init}; // valor inicial del fit (opcionalmente perturbado)
    if(jitter > 0)
    {
        for(auto& [key, v] : fitInit)
        {
            v[0] *= std::max(0.1, rng.Gaus(1, jitter));
            if(v.size() > 1 && key != "ps0")
                v[1] += rng.Gaus(0, 0.2 * jitter);
        }
    }
    const std::string outfile {"./Outputs/toy_fit_lDependent.root"};
    Fitters::RunFit(hToy, exmin, exmax, model, fitInit, inter.GetBounds(), inter.GetFixed(), outfile, "Toy fit",
                    {{"g0", "g.s"}, {"g1", "1st ex"}, {"v0", "3rd ex"}, {"ps0", "1-n phase"}}, false);

    // ---------------- Comprobaciones tras el fit ----------------
    // 1) El orden de ParNames coincide con el que he supuesto?
    {
        auto f {std::make_unique<TFile>(outfile.c_str())};
        auto* names {f->Get<std::vector<std::string>>("ParNames")};
        bool ok {names && names->size() == keyOfPar.size()};
        for(size_t i = 0; ok && i < names->size(); i++)
            ok = names->at(i).substr(0, names->at(i).find_first_of("_")) == keyOfPar[i];
        std::cout << (ok ? "Orden de parametros OK"
                         : "!! El orden de ParNames NO coincide con el supuesto: el toy esta mal generado")
                  << '\n';
    }
    // 2) Generado vs ajustado, por estado
    auto [vals, uncs] = Fitters::ReadInit(outfile);
    std::cout << std::setw(6) << "key" << std::setw(4) << "i" << std::setw(12) << "generated" << std::setw(12) << "fit"
              << std::setw(12) << "err" << std::setw(10) << "pull" << '\n';
    for(const auto& key : keys)
    {
        for(size_t i = 0; i < init.at(key).size(); i++)
        {
            double gen {init.at(key)[i]};
            double fit {vals.at(key)[i]};
            double err {uncs.at(key)[i]};
            std::cout << std::setw(6) << key << std::setw(4) << i << std::setw(12) << gen << std::setw(12) << fit
                      << std::setw(12) << err;
            if(err > 0)
                std::cout << std::setw(10) << std::setprecision(3) << (fit - gen) / err << std::setprecision(6);
            else
                std::cout << std::setw(10) << "fixed";
            std::cout << '\n';
        }
    }
}