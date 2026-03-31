#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActParticle.h"
#include "ActSRIM.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TH2D.h"
#include "TString.h"

#include <string>
#include <vector>

#include "../../PrettyStyle.C"


void pidFrontEachSilicon()
{
    PrettyStyle(true, true);

    std::string beam {"7Li"};
    std::string target {"d"};
    std::string light {"t"};
    std::string heavy {"7Li"};

    ActPhysics::Particle pb {beam};
    ActPhysics::Particle pt {target};
    ActPhysics::Particle pl {light};

    // --------------------------------------------------
    // Get the data for all runs to see pid for each silicon
    // --------------------------------------------------

    std::string dataconf {};
    if(beam == "11Li")
        dataconf = "../../configs/data_11Li.conf";
    else if(beam == "7Li")
        dataconf = "../../configs/data_7Li.conf";
    else
        throw std::runtime_error("Beam cannot differ from 11Li or 7Li");

    // Read data
    ActRoot::DataManager dataman {dataconf, ActRoot::ModeType::EMerge};
    auto chain {dataman.GetChain()};

    // RDataFrame
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame dforigin {*chain};

    // Filter silicon pads
    auto df = dforigin.Filter(
        [](ActRoot::MergerData& m)
        {
            // Mask L0_9
            if(m.fLight.fNs.size())
                if((m.fLight.fLayers.front() == "l0"))
                    return false;
            // Mask F0_2
            if(m.fLight.fNs.size())
                if((m.fLight.fLayers.front() == "f0") && (m.fLight.fNs.front() == 2))
                    return false;
            // Mask R0_3
            if(m.fLight.fNs.size())
                if((m.fLight.fLayers.front() == "r0"))
                    return false;

            if(m.fRun > 35 && m.fRun < 45)
            {
                if(!m.fLight.fLayers.empty() && m.fLight.fLayers.front() == "f0")
                {
                    if(!m.fLight.fNs.empty() && m.fLight.fNs.front() == 5)
                    {
                        return false;
                    }
                    else
                        return true;
                }
                else
                    return true;
            }
            else
                return true;
        },
        {"MergerData"});

    // Only f0 events
    auto def = df.Filter([](ActRoot::MergerData& m)
                         { return m.fLight.GetNLayers() == 1 && m.fLight.GetLayer(0) == "f0"; }, {"MergerData"});

    // ------------------------------------------------
    // Filter the df for each silicon pid
    // ------------------------------------------------
    // Strat reading the cuts for each silicon
    ActRoot::CutsManager<std::string> cuts;
    // Gas PID
    for(int i = 0; i < 12; i++)
    {
        cuts.ReadCut(TString::Format("f0_%d", i).Data(),
                     TString::Format("./Cuts/pid_%s_f0_%d_%s.root", light.c_str(), i, beam.c_str()).Data());
    }
    // Filter the df taking into account the hit in the silcion and whether it is inside the cut
    auto dfF0each = def.Filter(
        [&cuts](ActRoot::MergerData& m)
        {
            if(m.fLight.GetNLayers() != 1)
                return false;
            if(m.fLight.GetLayer(0) != "f0")
                return false;
            if(m.fLight.fNs.empty())
                return false;
            if(m.fLight.fEs.empty())
                return false;

            int sil = m.fLight.fNs.front();
            if(sil < 0 || sil >= 12)
                return false;

            std::string key {TString::Format("f0_%d", sil).Data()};
            if(!cuts.GetCut(key))
                return false;
            return cuts.IsInside(key, m.fLight.fEs.front(), m.fLight.fQave);
        },
        {"MergerData"});

    // Create EVErtex to plot the kinematic lines
    auto* srim {new ActPhysics::SRIM};
    // Correct SRIM names
    std::string srimName {};
    if(light == "p")
        srimName = "1H";
    if(light == "d")
        srimName = "2H";
    if(light == "t")
        srimName = "3H";
    srim->ReadTable(light, TString::Format("../../Calibrations/SRIM/%s_900mb_CF4_95-5.txt", srimName.c_str()).Data());
    srim->ReadTable(beam, TString::Format("../../Calibrations/SRIM/%s_900mb_CF4_95-5.txt", beam.c_str()).Data());
    srim->ReadTable("mylar", TString::Format("../../Calibrations/SRIM/%s_Mylar.txt", beam.c_str()).Data());
    // Build energy at vertex
    auto dfVertex = dfF0each.Define("EVertex",
                                    [&](const ActRoot::MergerData& d)
                                    {
                                        double ret {};
                                        if(d.fLight.IsFilled())
                                            ret = srim->EvalInitialEnergy(light, d.fLight.fEs.front(), d.fLight.fTL);
                                        else // L1 trigger
                                            ret = srim->EvalEnergy(light, d.fLight.fTL);
                                        return ret;
                                    },
                                    {"MergerData"});
    double initialEnergy {7.558}; // meassured by operators; resolution of 0,19%
    initialEnergy = srim->Slow("mylar", initialEnergy * pb.GetAMU(), 0.0168);
    initialEnergy = srim->Slow(beam, initialEnergy, 60); // 60 mm of gas before the pad plane
    initialEnergy = initialEnergy / pb.GetAMU();
    auto dfBeam = dfVertex.Define("EBeam", [&](const ActRoot::MergerData& d)
                                  { return srim->Slow(beam, initialEnergy * pb.GetAMU(), d.fRP.X()); }, {"MergerData"});

    ActPhysics::Kinematics kinTheo {pb, pt, pl, initialEnergy * pb.GetAMU()}; // energy meassured by operators
    ActPhysics::Kinematics kin {pb, pt, pl, initialEnergy * pb.GetAMU()};     // energy meassured by operators
    // Vector of kinematics as one object is needed per
    // processing slot (since we are changing EBeam in each entry)
    std::vector<ActPhysics::Kinematics> vkins {def.GetNSlots()};
    for(auto& k : vkins)
        k = kin;
    auto dfEx = dfBeam.DefineSlot("Ex",
                                  [&](unsigned int slot, const ActRoot::MergerData& d, double EVertex, double EBeam)
                                  {
                                      vkins[slot].SetBeamEnergy(EBeam);
                                      return vkins[slot].ReconstructExcitationEnergy(EVertex, (d.fThetaLight) *
                                                                                                  TMath::DegToRad());
                                  },
                                  {"MergerData", "EVertex", "EBeam"});

    // Gate in the heavy detection to ensure we are looking at the right reaction channel
    for(int i = 0; i < 4; i++)
    {
        cuts.ReadCut(
            TString::Format("f2_%d", i).Data(),
            TString::Format("../../PostAnalysis/Cuts/pid_%s_f2_%d_%s.root", heavy.c_str(), i, beam.c_str()).Data());
    }
    // Filter the df for the heavy
    auto dfExHeavy = dfEx.Filter(
        [&cuts](ActRoot::MergerData& m)
        {
            auto key {"f2_" + std::to_string(m.fHeavy.fNs[0])};
            if(cuts.IsInside(key, m.fHeavy.fEs[1], m.fHeavy.fEs[0]))
                return true;
            return false;
        },
        {"MergerData"});


    // --------------------------------------------------
    // Create histograms (one per silicon)
    // --------------------------------------------------
    std::vector<TH2D*> hSil(12, nullptr);

    for(int i = 0; i < 12; i++)
    {
        hSil[i] = new TH2D(TString::Format("hPID_sil_%d", i),
                           TString::Format("PID silicon %d;E_{Sil} [MeV];Q_{ave}", i), 450, 0, 55, 60, 0, 500);
    }

    std::vector<TH2D*> hSil_heavy(4, nullptr);
    for(int i = 0; i < 4; i++)
    {
        hSil_heavy[i] =
            new TH2D(TString::Format("hPID_heavy_f2_%d", i),
                     TString::Format("PID heavy f2 pad %d;E_{Heavy} [MeV];Q_{Heavy}", i), 800, 0, 80, 800, 0, 40);
    }

    std::vector<TH2D*> hSil_heavyCut(4, nullptr);
    for(int i = 0; i < 4; i++)
    {
        hSil_heavyCut[i] =
            new TH2D(TString::Format("hPIDcut_heavy_f2_%d", i),
                     TString::Format("PID heavy f2 pad %d cut in PID light;E_{Heavy} [MeV];Q_{Heavy}", i), 800, 0, 80,
                     800, 0, 40);
    }

    // Create histogram for kinematic lines
    ROOT::RDF::TH2DModel Kin_total {
        "hKin_total",
        TString::Format("Kinematics for %s ;#theta_{Lab} [#circ];E_{Vertex} [MeV]", light.c_str()),
        160,
        0,
        80,
        220,
        0,
        50};
    ROOT::RDF::TH1DModel Ex {"hEx",
                             TString::Format("Excitation energy for %s ;E_{x} [MeV];Counts / %.f keV", light.c_str(),
                                             (30. - (-15.)) / 300 * 1e3),
                             300, -15, 30};
    ROOT::RDF::TH2DModel ThetaLight_ThetaHeavy {
        "hThetaLight_ThetaHeavy",
        TString::Format("Theta light vs Theta heavy for %s;#theta_{Light} [#circ];#theta_{Heavy} [#circ]",
                        light.c_str()),
        120,
        0,
        60,
        120,
        0,
        30};


    // --------------------------------------------------
    // Fill histograms
    // --------------------------------------------------
    def.Foreach(
        [&hSil, &hSil_heavy](const ActRoot::MergerData& mer)
        {
            // Safety checks
            if(mer.fLight.fNs.empty())
                return;
            if(mer.fLight.fEs.empty())
                return;

            int sil = mer.fLight.fNs.front();
            if(sil < 0 || sil >= 12)
                return;

            hSil[sil]->Fill(mer.fLight.fEs.front(), mer.fLight.fQave);

            if(mer.fHeavy.GetNLayers() == 2)
            {
                int f2pad = mer.fHeavy.fNs[0];
                if(f2pad >= 0 && f2pad < 4)
                    hSil_heavy[f2pad]->Fill(mer.fHeavy.fEs[1], mer.fHeavy.fEs[0]);
            }
        },
        {"MergerData"});

    dfExHeavy.Foreach(
        [&hSil_heavyCut](const ActRoot::MergerData& mer)
        {
            if(mer.fHeavy.GetNLayers() == 2)
            {
                int f2pad = mer.fHeavy.fNs[0];
                if(f2pad >= 0 && f2pad < 4)
                    hSil_heavyCut[f2pad]->Fill(mer.fHeavy.fEs[1], mer.fHeavy.fEs[0]);
            }
        },
        {"MergerData"});

    auto hKin_total = dfEx.Histo2D(Kin_total, "fThetaLight", "EVertex");
    auto hEx = dfEx.Histo1D(Ex, "Ex");
    auto hThetaLight_ThetaHeavy = dfEx.Histo2D(ThetaLight_ThetaHeavy, "fThetaLight", "fThetaHeavy");
    // Get histos for heavy gated
    auto hKin_total_heavy = dfExHeavy.Histo2D(Kin_total, "fThetaLight", "EVertex");
    auto hEx_heavy = dfExHeavy.Histo1D(Ex, "Ex");
    auto hThetaLight_ThetaHeavy_heavy = dfExHeavy.Histo2D(ThetaLight_ThetaHeavy, "fThetaLight", "fThetaHeavy");


    // --------------------------------------------------
    // Draw
    // --------------------------------------------------
    auto* c = new TCanvas("c_pid_f0", "PID for f0", 1200, 800);
    c->Divide(4, 3);

    for(int i = 0; i < 12; i++)
    {
        c->cd(i + 1);
        hSil[i]->Draw("colz");
        if(cuts.GetCut(TString::Format("f0_%d", i).Data()))
            cuts.DrawCut(TString::Format("f0_%d", i).Data());
    }

    auto* cHeavyPID = new TCanvas("cHeavyPID", "Heavy PID for f0", 1200, 400);
    cHeavyPID->Divide(4, 2);
    for(int i = 0; i < 8; i++)
    {
        if(i < 4)
        {
            cHeavyPID->cd(i + 1);
            hSil_heavy[i]->Draw("colz");
            if(cuts.GetCut(TString::Format("f2_%d", i).Data()))
                cuts.DrawCut(TString::Format("f2_%d", i).Data());
        }
        else
        {
            cHeavyPID->cd(i + 1);
            hSil_heavyCut[i - 4]->Draw("colz");
            if(cuts.GetCut(TString::Format("f2_%d", i - 4).Data()))
                cuts.DrawCut(TString::Format("f2_%d", i - 4).Data());
        }
    }

    auto* cKin = new TCanvas("cKin", "Kinematics for f0", 800, 600);
    cKin->Divide(3, 1);
    cKin->cd(1);
    hKin_total->DrawClone("colz");
    kinTheo.GetKinematicLine3()->DrawClone("same");
    cKin->cd(2);
    hEx->DrawClone();
    cKin->cd(3);
    hThetaLight_ThetaHeavy->DrawClone("colz");
    kinTheo.GetTheta3vs4Line()->DrawClone("same");

    auto* cHeavy = new TCanvas("cHeavy", "Heavy gated for f0", 800, 600);
    cHeavy->Divide(3, 1);
    cHeavy->cd(1);
    hKin_total_heavy->DrawClone("colz");
    kinTheo.GetKinematicLine3()->DrawClone("same");
    cHeavy->cd(2);
    hEx_heavy->DrawClone();
    cHeavy->cd(3);
    hThetaLight_ThetaHeavy_heavy->DrawClone("colz");
    kinTheo.GetTheta3vs4Line()->DrawClone("same");
}
