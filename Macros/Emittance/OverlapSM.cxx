#include "ActSilMatrix.h"
#include "ActSilSpecs.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TPaveText.h"
#include "TString.h"

#include <iostream>
#include <memory>
#include <numeric>
#include <set>
#include <string>
#include <vector>

void OverlapSM()
{
    std::string beam {"7Li"};

    // Read histograms
    std::string filename {"./Outputs/histos" + beam + ".root"};
    auto file {new TFile {filename.c_str()}};
    auto* hxz {file->Get<TH2D>("hTrajXZ")};
    auto* hyz {file->Get<TH2D>("hTrajYZ")};

    // Mean of beam in histogram
    auto meanSide {hxz->GetMean(2)};
    auto meanFront {hyz->GetMean(2)};

    std::string label {"f0"};
    std::string layer {"f0"};

    // Real silicon specs
    ActPhysics::SilSpecs specs;
    specs.ReadFile("../../configs/silspecs.conf");

    // Read SM
    std::string filenameSM {"../SilVetos/Outputs/Dists/sms_" + label + ".root"};
    auto fileSM {new TFile {filenameSM.c_str()}};

    std::string name {};
    if(layer == "l0")
        name = "sm5";
    else if(layer == "r0")
        name = "sm6";
    else if(layer == "f0")
        name = "sm3";

    ActPhysics::SilMatrix* sm = fileSM->Get<ActPhysics::SilMatrix>(name.c_str());
    if(!sm)
    {
        std::cerr << "Error: could not find object " << name << " in file.\n";
        return;
    }

    // Clone to physical SM
    auto* phys = sm->Clone();
    phys->SetName(label);

    // Select silicons and reference depending on layer
    TH2D* h {};
    double ref {};
    std::set<int> sils;

    if(layer == "l0")
    {
        sils = {4, 5};
        h = hxz;
        ref = meanSide;
    }
    else if(layer == "r0")
    {
        sils = {6, 7};
        h = hxz;
        ref = meanSide;
    }
    else if(layer == "f0")
    {
        sils = {6, 7, 4};
        h = hyz;
        ref = meanFront;
    }

    // Compute mean Z of selected silicons
    std::vector<double> temp;
    for(auto sil : sils)
    {
        double x {}, y {};
        sm->GetSil(sil)->Center(x, y);
        temp.push_back(y);
    }

    double zmean = std::accumulate(temp.begin(), temp.end(), 0.0) / temp.size();
    double diff = zmean - ref;

    auto* text {new TPaveText {0.4, 0.75, 0.6, 0.88, "NDC"}};
    text->AddText(TString::Format("#DeltaZ = %.2f mm", diff));
    text->SetBorderSize(0);

    std::cout << "Mean for " << label << " : " << zmean << " mm" << '\n';

    // Move physical SM
    phys->MoveZTo(zmean, sils);

    // === Draw ===
    auto* c0 {new TCanvas {"c0", "SM and Emittance canvas"}};
    h->DrawCopy()->SetTitle(label.c_str());
    sm->Draw("same");
    text->Draw();

    c0->Modified();
    c0->Update();

    auto* c1 {new TCanvas {"c1", "Physical silicons"}};
    phys->Draw(false);
    auto* cl {sm->Clone()};
    cl->SetSyle(false, 0, 0, 3001);
    cl->Draw();
}
