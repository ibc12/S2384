#include "ActColors.h"
#include "ActKinematics.h"
#include "ActParticle.h"

#include "ROOT/RDataFrame.hxx"
#include "Rtypes.h"

#include "TCanvas.h"
#include "TEfficiency.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"
#include "TVirtualPad.h"

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../PostAnalysis/HistConfig.h"

void Plotter(const std::vector<double>& Exs, const std::string& beam, const std::string& target,
             const std::string& light, double T1, int neutronPS, int protonPS)
{
    ROOT::EnableImplicitMT();

    // Set if we should fit or not
    bool isPS {(neutronPS != 0 || protonPS != 0)};
    if(!isPS)
        std::cout << BOLDCYAN << "Plotting usual simulation" << RESET << '\n';
    else
        std::cout << BOLDCYAN << "Plotting a phase space" << RESET << '\n';

    // Save histograms
    std::vector<TH1D*> hsEx, hsRPx, hsCM;
    std::vector<TH2D*> hsKin, hsSP, hsRP, hsRPCM;
    std::vector<TEfficiency*> effsCM, effsLab;
    // Iterate
    int idx {1};
    for(const auto& Ex : Exs)
    {
        TString fileName {TString::Format("./Outputs/%s/%s_%s_TRIUMF_Eex_%.3f_nPS_%d_pPS_%d.root", beam.c_str(),
                                          target.c_str(), light.c_str(), Ex, neutronPS, protonPS)};
        // Get DF (this is fine to construct while MT is enabled)
        ROOT::RDataFrame df("SimulationTTree", fileName);

        // Book histograms (these evaluate on the dataframe and return owned histograms)
        auto hEx {
            df.Histo1D(HistConfig::ChangeTitle(HistConfig::Ex, TString::Format("%s(%s, %s) Ex = %.2f", beam.c_str(),
                                                                               target.c_str(), light.c_str(), Ex)),
                       "Eex", "weight")};
        auto hKin {df.Histo2D(
            HistConfig::ChangeTitle(HistConfig::KinSimu, TString::Format("%s(%s, %s) Ex = %.2f", beam.c_str(),
                                                                         target.c_str(), light.c_str(), Ex)),
            "theta3Lab", "EVertex", "weight")};
        auto hCM {df.Histo1D(HistConfig::ThetaCM, "theta3CM")};
        hCM->Scale(1. / (hCM->Integral() > 0 ? hCM->Integral() : 1.0));
        hCM->SetTitle(TString::Format("%.2f MeV;#theta_{CM} [#circ];Normalized counts", Ex));

        // Read detectors/efficiencies/histos directly from file
        auto* f {new TFile(fileName, "READ")};

        // --- Clone TEfficiency objects so they survive file close ---
        auto* effCM_raw = f->Get<TEfficiency>("effCM");
        if(!effCM_raw) {
            f->Close();
            delete f;
            throw std::runtime_error(std::string("Could not read TEfficiency named effCM in file ") + std::string(fileName.Data()));
        }
        TEfficiency* effCM_clone = dynamic_cast<TEfficiency*>(effCM_raw->Clone());
        effCM_clone->SetDirectory(nullptr); // detach from any file

        auto* effLab_raw = f->Get<TEfficiency>("effLab");
        if(!effLab_raw) {
            f->Close();
            delete f;
            throw std::runtime_error(std::string("Could not read TEfficiency named effLab in file ") + std::string(fileName.Data()));
        }
        TEfficiency* effLab_clone = dynamic_cast<TEfficiency*>(effLab_raw->Clone());
        effLab_clone->SetDirectory(nullptr);

        // Get hRP (we keep it as before)
        auto* hRP {f->Get<TH2D>("hRP")};
        if(!hRP) {
            f->Close();
            delete f;
            throw std::runtime_error(std::string("Could not read TH2D named hRP in file ") + std::string(fileName.Data()));
        }
        hRP->SetDirectory(nullptr);

        // Close file now that we've cloned/detached everything we need
        f->Close();
        delete f;

        // clone in order to save
        hsEx.push_back((TH1D*)hEx->Clone());
        hsKin.push_back((TH2D*)hKin->Clone());
        hsCM.push_back((TH1D*)hCM->Clone());
        effsCM.push_back(effCM_clone);
        effsLab.push_back(effLab_clone);
        hsRP.push_back(hRP);
        idx++;
    }
    // Fit to gaussians!
    auto* gsigmas {new TGraphErrors()};
    double range {4.5}; // MeV around the mean
    for(int i = 0; i < (int)hsEx.size(); i++)
    {
        if(isPS)
            continue;
        hsEx[i]->Fit("gaus", "0Q", "", Exs[i] - range, Exs[i] + range);
        auto* f {hsEx[i]->GetFunction("gaus")};
        if(!f)
            continue;
        gsigmas->SetPoint(gsigmas->GetN(), Exs[i], f->GetParameter("Sigma"));
        gsigmas->SetPointError(gsigmas->GetN() - 1, 0, f->GetParError(2));
    }

    // --- Disable implicit MT before any drawing / canvas creation (graphics are not thread-safe) ---
    ROOT::DisableImplicitMT();

    // Plot!
    std::vector<TCanvas*> cs;
    for(int i = 0; i < (int)hsEx.size(); i++)
    {
        cs.push_back(new TCanvas {TString::Format("c%d", i), TString::Format("Ex = %.2f MeV", Exs[i])});
        cs[i]->DivideSquare(6);
        // Get theoretical kinematics
        ActPhysics::Particle p1 {beam};
        double T1Total {T1};
        ActPhysics::Kinematics kin(beam, target, light, T1Total, Exs[i]);
        auto* gtheo {kin.GetKinematicLine3()};
        // Kinematics
        cs[i]->cd(1);
        hsKin[i]->Draw("colz");
        if(gtheo) gtheo->Draw("same");
        // Ex
        cs[i]->cd(2);
        if(!isPS)
            hsEx[i]->GetXaxis()->SetRangeUser(Exs[i] - range, Exs[i] + range);
        hsEx[i]->Draw("hist");
        if (hsEx[i]->GetListOfFunctions()) {
            for(auto* o : *(hsEx[i]->GetListOfFunctions()))
                if(o)
                    o->Draw("same");
        }
        // Get a line centered at Ex
        gPad->Update();
        auto* line {new TLine(Exs[i], gPad->GetUymin(), Exs[i], gPad->GetUymax())};
        line->SetLineWidth(2);
        line->SetLineColor(kMagenta);
        line->Draw();
        // Efficiency
        cs[i]->cd(3);
        if(effsCM[i]) effsCM[i]->Draw("apl");
        // RPx distribution -> using hsRP (you previously referenced hsRPx but never filled it)
        cs[i]->cd(4);
        // If you want a 1D projection for "RPx", create it now (example projection on X)
        TH1* projX = hsRP[i]->ProjectionX(); // owns a new histogram
        if(projX) { projX->Draw("hist"); }
        cs[i]->cd(5);
        hsRP[i]->Draw("colz"); // instead of hsRPCM which wasn't filled
        cs[i]->cd(6);
        if(effsLab[i]) effsLab[i]->Draw("apl");
        // hsCM[i]->Draw("hist");
        // hCMExp->Draw("hist same");
    }

    if(isPS)
        return;

    auto* csigma {new TCanvas("csigma", "Sigmas from fits")};
    csigma->DivideSquare(1 + hsRP.size());
    csigma->cd(1);
    gsigmas->SetTitle(";E_{x} [MeV];#sigma in E_{x} [MeV]");
    gsigmas->SetMarkerStyle(24);
    gsigmas->SetLineColor(kViolet);
    gsigmas->SetLineWidth(2);
    gsigmas->Draw("apl0");
    for(int i = 0; i < (int)hsRP.size(); i++)
    {
        csigma->cd(2 + i);
        // draw projection or histogram: reuse projX if you stored it, here we project again
        TH1* projX = hsRP[i]->ProjectionX();
        if(projX) projX->Draw("hist");
    }

    auto* file {new TFile {TString::Format("/media/Data/E796v2/Simulation/Outputs/%s/sigmas_%s_%s_%s.root",
                                           beam.c_str(), beam.c_str(), target.c_str(), light.c_str()),
                           "recreate"}};
    file->cd();
    gsigmas->Write("gsigmas");
    file->Close();
    delete file;
}
