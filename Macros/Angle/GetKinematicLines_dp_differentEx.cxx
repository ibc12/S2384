#include "ActDataManager.h"
#include "ActKinematicGenerator.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActTPCData.h"

#include <ROOT/RDataFrame.hxx>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TTree.h"
#include "TLine.h"
#include "TMultiGraph.h"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

void GetKinematicLines_dp_differentEx()
{
    std::string beam {"7Li"};
    std::string target {"2H"};
    std::string light {"1H"};
    std::string heavy {"8Li"};
    double Tbeam = 7 * 7.558; // MeV
    std::vector<double> Exs {0, 0.981, 2.255, 3.210, 4.100, 5.400, 6.100, 6.530, 7.100, 8.000, 9.000, 10.000}; // MeV

    // Colores distinguibles (evito amarillo, que se ve mal sobre fondo blanco)
    std::vector<int> colors {kBlack,     kRed + 1,  kOrange + 7, kGreen + 2,  kTeal + 3,
                             kAzure + 1, kBlue + 1, kViolet + 1, kMagenta + 1};

    std::vector<ActPhysics::Kinematics*> kinVec;
    for(const auto& Ex : Exs)
        kinVec.push_back(new ActPhysics::Kinematics {beam, target, light, heavy, Tbeam, Ex});

    auto* canvas = new TCanvas("cKinematicLines", "Kinematic lines for 7Li(d,p) at different Ex", 1400, 1000);
    auto* mg = new TMultiGraph();
    mg->SetTitle("^{7}Li(d,p) kinematic lines;#theta_{Lab} [#circ];E_{p} [MeV]");
    auto* leg = new TLegend(0.70, 0.55, 0.89, 0.89);
    leg->SetHeader("E_{x} [MeV]", "C");
    leg->SetBorderSize(0);

    for(size_t i = 0; i < kinVec.size(); ++i)
    {
        auto* g = kinVec[i]->GetKinematicLine3();
        g->SetLineColor(colors[i % colors.size()]);
        g->SetLineWidth(2);
        mg->Add(g, "L");
        leg->AddEntry(g, TString::Format("%.3f", Exs[i]), "l");
    }
    mg->Draw("A");
    leg->Draw();

    // Opcional: ventana de energía de tu detector (0.5 - 2 MeV)
    for(double e : {0.5, 2.0})
    {
        auto* l = new TLine(mg->GetXaxis()->GetXmin(), e, mg->GetXaxis()->GetXmax(), e);
        l->SetLineStyle(2);
        l->SetLineColor(kGray + 2);
        l->Draw();
    }
}