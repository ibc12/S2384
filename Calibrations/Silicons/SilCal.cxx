#include "ActSRIM.h"

#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TVirtualPad.h"

#include "CalibrationRunner.h"
#include "CalibrationSource.h"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

const std::string layer {"f0"};

std::vector<TH1D*> ReadData(const std::string& file, const std::string& dir, const std::string& label)
{
    auto* f {new TFile {file.c_str()}};
    // f->ls();
    auto updir {f->Get<TDirectory>("Raw")};
    updir->ls();
    auto lowdir {updir->Get<TDirectory>(dir.c_str())};
    // lowdir->ls();
    // Read
    std::vector<TH1D*> ret;
    auto keys {lowdir->GetListOfKeys()};
    for(auto* key : *keys)
    {
        std::string str {key->GetName()};
        auto idx {str.find_first_of("0123456789")};
        auto name {str.substr(0, idx)};
        if(!(name == label))
            continue;
        ret.push_back((TH1D*)lowdir->Get<TH1I>(key->GetName()));
    }
    return ret;
}

void CorrectSource(Calibration::Source* source, ActPhysics::SRIM* srim, const std::string& table, double thickness,
                   double angle = 0)
{
    auto& energies {source->GetRefToEnergies()};
    auto& limits {source->GetRefToLimits()};
    auto labels {source->GetLabels()};
    // Peak energies
    for(auto& source : energies)
        for(auto& energy : source)
            energy = srim->Slow(table, energy, thickness, angle);
    // Peak boundaries
    for(auto& [key, vals] : limits)
    {
        vals.first = srim->Slow(table, vals.first, thickness, angle);
        vals.second = srim->Slow(table, vals.second, thickness, angle);
    }
}

void e864()
{
    std::string which {"fu"};
    std::string label {};
    if(which == "fu")
        label = "SI_";
    else if(which == "fc")
        label = "SI_MGP_";
    else
        throw std::runtime_error("Which can only be fu or fc");
    // Read data
    // auto hs {ReadData("./Inputs/e864/Si_USC_for_e864_afterBOOM.root", "SI", label)};
    auto hs {ReadData("./Inputs/e864/Si_USC_for_e864_6_6_24.root", "SI", label)};
    // Pick only necessary
    // hs = {hs[0], hs[1]};

    // Source of ganil
    Calibration::Source source {};
    // source.Print();

    // Correct by energy losses in Al dead layer
    ActPhysics::SRIM srim {"al", "./Inputs/alpha_Al.txt"};
    CorrectSource(&source, &srim, "al", 0.5e-3); // 0.5 um to mm
    source.Print();

    // Rebin
    std::vector<TH1D*> hsrebin;
    int idx {};
    for(auto& h : hs)
    {
        h->Rebin(16);
        hsrebin.push_back((TH1D*)h->Clone());
        idx++;
    }
    // Runner per silicon
    std::vector<Calibration::Runner> runners;
    // Graph
    auto* gr {new TGraphErrors};
    gr->SetNameTitle("g", "Resolution;Silicon index;#sigma ^{241}Am [keV]");
    // Save
    std::ofstream streamer {"./Outputs/junk_e864_" + which + ".dat"};
    streamer << std::fixed << std::setprecision(8);
    std::vector<std::shared_ptr<TH1D>> hfs;
    for(int s = 0; s < hsrebin.size(); s++)
    {
        runners.emplace_back(&source, hsrebin[s], hs[s]);
        runners.back().SetGaussPreWidth(150);
        if(which == "fc")
            runners.back().SetRange(4500, 6500);
        else
            runners.back().SetRange(6000, 9500);
        // if(s == 1)
        runners.back().DisableXErrors();
        runners.back().DoIt();
        auto* c {new TCanvas};
        runners.back().Draw(c);
        std::cout << "Sil index : " << s << " hist name : " << hs[s]->GetName() << '\n';
        runners.back().PrintRes();
        std::cout << '\n';
        gr->SetPoint(s, s + 1, runners.back().GetRes("241Am") * 1e3);
        gr->SetPointError(s, 0, runners.back().GetURes("241Am") * 1e3);
        hfs.push_back(runners.back().GetHistFinal());

        // Save calibration in file
        auto label {TString::Format("Sil_%s_%d_E", which.c_str(), s)};
        auto labelp {TString::Format("Sil_%s_%d_P", which.c_str(), s)};
        auto [p0, p1] {runners.back().GetParameters()};
        streamer << label << " " << p0 << " " << p1 << '\n';
        auto [ped, peds] {runners.back().GetPedestal()};
        streamer << labelp << " " << ped << " " << peds << '\n';
    }
    streamer.close();

    // Plot
    auto* c0 {new TCanvas {"c0", "Raw silicon data"}};
    c0->DivideSquare(hs.size());
    for(int i = 0; i < hs.size(); i++)
    {
        c0->cd(i + 1);
        gPad->SetLogy();
        hs[i]->Draw();
    }

    auto* c1 {new TCanvas {"c11", "Resolution canvas"}};
    gr->SetMarkerStyle(25);
    gr->SetLineWidth(2);
    gr->Draw("apl");

    auto* c2 {new TCanvas {"c2", "Final his canvas"}};
    // gROOT->SetSelectedPad(nullptr);
    c2->DivideSquare(hfs.size());
    for(int p = 0; p < hfs.size(); p++)
    {
        c2->cd(p + 1);
        hfs[p]->SetTitle(TString::Format("%s%d", label.c_str(), p + 1));
        hfs[p]->DrawClone();
        for(auto* o : *(hfs[p]->GetListOfFunctions()))
            if(o)
                o->DrawClone("same");
    }
}
