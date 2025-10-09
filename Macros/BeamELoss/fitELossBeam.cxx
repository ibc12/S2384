#include "TFile.h"
#include "TH1.h"
#include "TH2D.h"
#include "TString.h"
#include "TCanvas.h"
#include "ActSRIM.h"
#include "TF1.h"
#include "ActLine.h"

void fitELossBeam()
{
    // Get histos from files
    TFile* f7Li = TFile::Open("./Outputs/DE_7Li_10pads.root");
    TFile* f11Li = TFile::Open("./Outputs/DE_11Li_10pads.root");
    TH1D* h7Li = f7Li->Get<TH1D>("hDE_7Li");
    TH1D* h11Li = f11Li->Get<TH1D>("hDE_11Li");

    auto srim {new ActPhysics::SRIM};
    srim->ReadTable("7Li", "../../Calibrations/SRIM/7Li_900mb_CF4_95-5.txt");
    srim->ReadTable("11Li", "../../Calibrations/SRIM/11Li_900mb_CF4_95-5.txt");
    double E_7Li {7.5 * 7};
    double E_11Li {7.5 * 11};

    auto eLoss_7Li {E_7Li - srim->Slow("7Li", E_7Li, 20)};
    auto eLoss_11Li {E_11Li - srim->Slow("11Li", E_11Li, 20)};

    // Fit the histograms into gausians
    h7Li->Fit("gaus", "Q");
    h11Li->Fit("gaus", "Q");
    // Get functions to plot them
    auto f_7Li {h7Li->GetFunction("gaus")};
    auto f_11Li {h11Li->GetFunction("gaus")};

    // Let's do two lines, starting in 0,0 and going to eLoss from srim, and the center of the gausian fit to calibrate
    auto mean_7Li {f_7Li->GetParameter(1)};
    auto mean_11Li {f_11Li->GetParameter(1)};

    ROOT::Math::XYZPointF point_0 {0,0,0};
    ROOT::Math::XYZPointF point_7Li {eLoss_7Li, mean_7Li, 0};
    ROOT::Math::XYZPointF point_11Li {eLoss_11Li, mean_11Li, 0};
    auto line_7Li {new ActRoot::Line(point_0, point_7Li)};
    auto line_11Li {new ActRoot::Line(point_0, point_11Li)};

    // Get the sigmas and convert it to energy
    auto sigma_7Li_raw {f_7Li->GetParameter(2)};
    auto sigma_11Li_raw {f_11Li->GetParameter(2)};
    auto sigma_7Li_MeV {line_7Li->MoveToY(sigma_7Li_raw).X()};
    auto sigma_11Li_MeV {line_11Li->MoveToY(sigma_11Li_raw).X()};

    // Do the couts for the mean and the sigmas
    std::cout << "7Li: E_loss from SRIM: " << eLoss_7Li << " MeV, mean from fit: " << mean_7Li << " u.a., sigma from fit: " << sigma_7Li_raw << " u.a. = " << sigma_7Li_MeV << " MeV" << std::endl;
    std::cout << "11Li: E_loss from SRIM: " << eLoss_11Li << " MeV, mean from fit: " << mean_11Li << " u.a., sigma from fit: " << sigma_11Li_raw << " u.a. = " << sigma_11Li_MeV << " MeV" << std::endl;
    
    // Plot
    auto c{new TCanvas("c")};
    c->DivideSquare(2);
    c->cd(1);
    h7Li->DrawClone();
    f_7Li->SetLineColor(kRed);
    f_7Li->Draw("same");
    c->cd(2);
    h11Li->DrawClone();
    f_11Li->SetLineColor(kRed);
    f_11Li->Draw("same");
}