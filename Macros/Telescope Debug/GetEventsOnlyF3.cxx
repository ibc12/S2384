#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActParticle.h"
#include "ActRunner.h"
#include "ActSRIM.h"
#include "ActSilData.h"
#include "ActSilMatrix.h"
#include "ActSilSpecs.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TH2.h"
#include "TH2D.h"
#include "TMath.h"
#include "TRandom.h"
#include "TString.h"

#include "Math/Point3Dfwd.h"

#include <fstream>
#include <map>
#include <string>

#include "../../PostAnalysis/HistConfig.h"


void GetEventsOnlyF3()
{
    std::string beam {"11Li"};
    // ROOT::EnableImplicitMT();
    ROOT::RDataFrame dfini {"PreProcessed_Tree",
                            TString::Format("../../PostAnalysis/Outputs/tree_preprocess_F_%s.root", beam.c_str())};
    auto df {dfini.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == true; }, {"MergerData"})};

    // histograms for hit in f3 but no f2
    auto hSPf2_f3_no_f2_postL0 {HistConfig::SP.GetHistogram()};
    hSPf2_f3_no_f2_postL0->SetTitle("SP for f2 when no hit f2 but f3 yes - post L0 change");
    auto hSPf2_f3_no_f2_preL0 {HistConfig::SP.GetHistogram()};
    hSPf2_f3_no_f2_preL0->SetTitle("SP for f2 when no hit f2 but f3 yes - pre L0 change");
    auto hSPf3_f3_no_f2_postL0 {HistConfig::SP.GetHistogram()};
    hSPf3_f3_no_f2_postL0->SetTitle("SP for f3 when no hit f2 but f3 yes - post L0 change");
    auto hSPf3_f3_no_f2_preL0 {HistConfig::SP.GetHistogram()};
    hSPf3_f3_no_f2_preL0->SetTitle("SP for f3 when no hit f2 but f3 yes - pre L0 change");
    auto hEnergy_f3_no_f2_postL0 {new TH1D("hEnergy_f3_no_f2_postL0",
                                           "Energy for f3 when no hit f2 but f3 yes - post L0 change;E_{Vertex} [MeV]",
                                           140, 0, 70)};
    auto hEnergy_f3_no_f2_preL0 {new TH1D("hEnergy_f3_no_f2_preL0",
                                          "Energy for f3 when no hit f2 but f3 yes - pre L0 change;E_{Vertex} [MeV]",
                                          140, 0, 70)};
    // histograms for hit in f2 but no f3
    auto hSPf2_f2_no_f3_postL0 {HistConfig::SP.GetHistogram()};
    hSPf2_f2_no_f3_postL0->SetTitle("SP for f2 when hit f2 but no hit f3 - post L0 change");
    auto hSPf2_f2_no_f3_preL0 {HistConfig::SP.GetHistogram()};
    hSPf2_f2_no_f3_preL0->SetTitle("SP for f2 when hit f2 but no hit f3 - pre L0 change");
    auto hSPf3_f2_no_f3_postL0 {HistConfig::SP.GetHistogram()};
    hSPf3_f2_no_f3_postL0->SetTitle("SP for f3 when hit f2 but no hit f3 - post L0 change");
    auto hSPf3_f2_no_f3_preL0 {HistConfig::SP.GetHistogram()};
    hSPf3_f2_no_f3_preL0->SetTitle("SP for f3 when hit f2 but no hit f3 - pre L0 change");
    auto hEnergy_f2_no_f3_postL0 {new TH1D("hEnergy_f2_no_f3_postL0",
                                           "Energy for f2 when hit f2 but no hit f3 - post L0 change;E_{Vertex} [MeV]",
                                           100, 0, 20)};
    auto hEnergy_f2_no_f3_preL0 {new TH1D("hEnergy_f2_no_f3_preL0",
                                          "Energy for f2 when hit f2 but no hit f3 - pre L0 change;E_{Vertex} [MeV]",
                                          100, 0, 20)};
    auto hSize_f2Hits_postL0 {new TH1D("hSize_f2Hits_postL0", "Number of f2 hits - post L0 change;N_{f2 hits}", 4, 0, 3)};
    auto hSize_f2Hits_preL0 {new TH1D("hSize_f2Hits_preL0", "Number of f2 hits - pre L0 change;N_{f2 hits}", 4, 0, 3)};

    // Get silCentre points for post and pre
    std::string filenameSMfront_post {"../SilVetos/Outputs/Dists/sms_f0.root"};
    auto fileSMfront_post {new TFile {filenameSMfront_post.c_str()}};
    ActPhysics::SilMatrix* smfront_post =
        fileSMfront_post->Get<ActPhysics::SilMatrix>("sm7"); // matrix for good distance of left wall
    double silCentreFront_post = smfront_post->GetMeanZ({6, 7, 4});
    // Now pre
    std::string filenameSMfront_pre {"../SilVetos/Outputs/Dists/sms_f0_preL0change.root"};
    auto fileSMfront_pre {new TFile {filenameSMfront_pre.c_str()}};
    ActPhysics::SilMatrix* smfront_pre =
        fileSMfront_pre->Get<ActPhysics::SilMatrix>("sm5"); // matrix for good distance of left wall
    double silCentreFront_pre = smfront_pre->GetMeanZ({6, 7, 4});

    // Silspecs for post
    auto* sils {new ActPhysics::SilSpecs};
    std::string silConfig("silspecs"); // no front silicons, only lateral ones
    sils->ReadFile("../../configs/" + silConfig + ".conf");
    // double silCentreFront = 266.318;
    for(auto& [name, layer] : sils->GetLayers())
    {
        // Move the layers to the correct positions
        if(name == "f2")
            layer.MoveZTo(silCentreFront_post + 12.5,
                          {0}); // Manually substract 12.5 to get the f2 in the same height as the f3
        if(name == "f3")
            layer.MoveZTo(silCentreFront_post, {0});
    }
    // Silspecs for pre
    auto* sils_pre {new ActPhysics::SilSpecs};
    std::string silConfig_pre("silspecs"); // no front silicons, only lateral ones
    sils_pre->ReadFile("../../configs/" + silConfig_pre + ".conf");
    // double silCentreFront = 266.318;
    for(auto& [name, layer] : sils_pre->GetLayers())
    {
        // Move the layers to the correct positions
        if(name == "f2")
            layer.MoveZTo(silCentreFront_pre + 12.5,
                          {0}); // Manually substract 12.5 to get the f2 in the same height as the f3
        if(name == "f3")
            layer.MoveZTo(silCentreFront_pre, {0});
    }
    std::cout << "Silicon front centre post L0 change at Z = " << silCentreFront_post << " mm" << std::endl;
    std::cout << "Silicon front centre pre L0 change at Z = " << silCentreFront_pre << " mm" << std::endl;

    // Just to get the heavy direction
    ActSim::Runner runner(nullptr, nullptr, gRandom, 0);

    // Filter for hit in layer f3 and no hit in layer f2
    auto dfF3 = df.Filter([](ActRoot::SilData& sil)
                          { return sil.fSiE["f2"].size() == 0 && sil.fSiE["f3"].size() != 0; }, {"SilData"});
    // Filter for fit in layer f2 and not in f3
    auto dfF2 = df.Filter([](ActRoot::SilData& sil)
                          { return sil.fSiE["f2"].size() != 0 && sil.fSiE["f3"].size() == 0; }, {"SilData"});
    // Filter for hit in both layers f2 and f3
    auto dfF2F3 = df.Filter([](ActRoot::SilData& sil)
                            { return sil.fSiE["f2"].size() != 0 && sil.fSiE["f3"].size() != 0; }, {"SilData"});

    std::cout << "Number of events with hit in f3 and no hit in f2: " << dfF3.Count().GetValue() << std::endl;
    std::cout << "Number of events with hit in f2 and no hit in f3: " << dfF2.Count().GetValue() << std::endl;
    std::cout << "Number of events with hit in both f2 and f3: " << dfF2F3.Count().GetValue() << std::endl;

    // And for all events fill theextrapolated to f3
    dfF3.Foreach(
        [&](ActRoot::SilData& sil, ActRoot::MergerData& m, ActRoot::TPCData& tpc)
        {
            int heavyIdx = m.fHeavyIdx;
            auto heavyDirFloat = tpc.fClusters[heavyIdx].GetLine().GetDirection().Unit();
            ROOT::Math::XYZVector heavyDir {heavyDirFloat.X(), heavyDirFloat.Y(), heavyDirFloat.Z()};
            ActRoot::Line heavyLine {
                ROOT::Math::XYZPointF {static_cast<float>(m.fRP.X()), static_cast<float>(m.fRP.Y()),
                                       static_cast<float>(m.fRP.Z())},
                ROOT::Math::XYZVectorF {static_cast<float>(heavyDir.X()), static_cast<float>(heavyDir.Y()),
                                        static_cast<float>(heavyDir.Z())},
                -1};
            ROOT::Math::XYZPoint extrapolatedPointF3 {};
            ROOT::Math::XYZPoint extrapolatedPointF2 {};
            double energyF3 {};
            energyF3 = sil.fSiE["f3"].front();
            if(m.fRun > 66)
            {
                extrapolatedPointF3 = heavyLine.MoveToX(sils->GetLayer("f3").GetPoint().X());
                extrapolatedPointF2 = heavyLine.MoveToX(sils->GetLayer("f2").GetPoint().X());
                hSPf2_f3_no_f2_postL0->Fill(extrapolatedPointF2.Y(), extrapolatedPointF2.Z());
                hSPf3_f3_no_f2_postL0->Fill(extrapolatedPointF3.Y(), extrapolatedPointF3.Z());
                hEnergy_f3_no_f2_postL0->Fill(energyF3);
            }
            else
            {
                extrapolatedPointF3 = heavyLine.MoveToX(sils_pre->GetLayer("f3").GetPoint().X());
                extrapolatedPointF2 = heavyLine.MoveToX(sils_pre->GetLayer("f2").GetPoint().X());
                hSPf2_f3_no_f2_preL0->Fill(extrapolatedPointF2.Y(), extrapolatedPointF2.Z());
                hSPf3_f3_no_f2_preL0->Fill(extrapolatedPointF3.Y(), extrapolatedPointF3.Z());
                hEnergy_f3_no_f2_preL0->Fill(energyF3);
            }
        },
        {"SilData", "MergerData", "TPCData"});

    // Do the same for the case yes f2 but no f3
    dfF2.Foreach(
        [&](ActRoot::SilData& sil, ActRoot::MergerData& m, ActRoot::TPCData& tpc)
        {
            int heavyIdx = m.fHeavyIdx;
            auto heavyDirFloat = tpc.fClusters[heavyIdx].GetLine().GetDirection().Unit();
            ROOT::Math::XYZVector heavyDir {heavyDirFloat.X(), heavyDirFloat.Y(), heavyDirFloat.Z()};
            ActRoot::Line heavyLine {
                ROOT::Math::XYZPointF {static_cast<float>(m.fRP.X()), static_cast<float>(m.fRP.Y()),
                                       static_cast<float>(m.fRP.Z())},
                ROOT::Math::XYZVectorF {static_cast<float>(heavyDir.X()), static_cast<float>(heavyDir.Y()),
                                        static_cast<float>(heavyDir.Z())},
                -1};
            double energyF2 {};
            energyF2 = sil.fSiE["f2"].front();
            ROOT::Math::XYZPoint extrapolatedPointF3 {};
            ROOT::Math::XYZPoint extrapolatedPointF2 {};
            if(m.fRun > 66)
            {
                extrapolatedPointF3 = heavyLine.MoveToX(sils->GetLayer("f3").GetPoint().X());
                extrapolatedPointF2 = heavyLine.MoveToX(sils->GetLayer("f2").GetPoint().X());
                hSPf2_f2_no_f3_postL0->Fill(extrapolatedPointF2.Y(), extrapolatedPointF2.Z());
                hSPf3_f2_no_f3_postL0->Fill(extrapolatedPointF3.Y(), extrapolatedPointF3.Z());
                hEnergy_f2_no_f3_postL0->Fill(energyF2);
                hSize_f2Hits_postL0->Fill(sil.fSiE["f2"].size());
            }
            else
            {
                extrapolatedPointF3 = heavyLine.MoveToX(sils_pre->GetLayer("f3").GetPoint().X());
                extrapolatedPointF2 = heavyLine.MoveToX(sils_pre->GetLayer("f2").GetPoint().X());
                hSPf2_f2_no_f3_preL0->Fill(extrapolatedPointF2.Y(), extrapolatedPointF2.Z());
                hSPf3_f2_no_f3_preL0->Fill(extrapolatedPointF3.Y(), extrapolatedPointF3.Z());
                hEnergy_f2_no_f3_preL0->Fill(energyF2);
                hSize_f2Hits_preL0->Fill(sil.fSiE["f2"].size());
            }
        },
        {"SilData", "MergerData", "TPCData"});

    // Save events in a .dat file
    // std::ofstream outfileF3("./Outputs/events_f3_only.dat");
    // dfF3.Foreach([&outfileF3](ActRoot::MergerData m) { m.Stream(outfileF3); }, {"MergerData"});

    // Plot
    auto* cSP_f3_no_f2 {new TCanvas {"cSP_f3_no_f2", "Sil Points for events hit with f3 but no f2", 1200, 800}};
    cSP_f3_no_f2->DivideSquare(4);
    cSP_f3_no_f2->cd(1);
    hSPf2_f3_no_f2_postL0->DrawClone("colz");
    auto f2 = sils->GetLayer("f2").GetSilMatrix();
    f2->DrawClone();
    cSP_f3_no_f2->cd(2);
    hSPf3_f3_no_f2_postL0->DrawClone("colz");
    auto f3 = sils->GetLayer("f3").GetSilMatrix();
    f3->DrawClone();
    cSP_f3_no_f2->cd(3);
    hSPf2_f3_no_f2_preL0->DrawClone("colz");
    auto f2_pre = sils_pre->GetLayer("f2").GetSilMatrix();
    f2_pre->DrawClone();
    cSP_f3_no_f2->cd(4);
    hSPf3_f3_no_f2_preL0->DrawClone("colz");
    auto f3_pre = sils_pre->GetLayer("f3").GetSilMatrix();
    f3_pre->DrawClone();


    auto* cSP_f2_no_f3 {new TCanvas {"cSP_f2_no_f3", "Sil Points for events with f2 but no f3", 1200, 800}};
    cSP_f2_no_f3->DivideSquare(4);
    cSP_f2_no_f3->cd(1);
    hSPf2_f2_no_f3_postL0->DrawClone("colz");
    f2->DrawClone();
    cSP_f2_no_f3->cd(2);
    hSPf3_f2_no_f3_postL0->DrawClone("colz");
    f3->DrawClone();
    cSP_f2_no_f3->cd(3);
    hSPf2_f2_no_f3_preL0->DrawClone("colz");
    f2_pre->DrawClone();
    cSP_f2_no_f3->cd(4);
    hSPf3_f2_no_f3_preL0->DrawClone("colz");
    f3_pre->DrawClone();

    auto* cEnergy {new TCanvas {"cEnergy", "Energy for f2 and f3 when no hit in the other layer", 1200, 800}};
    cEnergy->DivideSquare(5);
    cEnergy->cd(1);
    hEnergy_f3_no_f2_preL0->DrawClone("colz");  
    cEnergy->cd(2);
    hEnergy_f3_no_f2_postL0->DrawClone("colz");
    cEnergy->cd(3);
    hEnergy_f2_no_f3_preL0->DrawClone("colz");
    cEnergy->cd(4);
    hEnergy_f2_no_f3_postL0->DrawClone("colz");
    cEnergy->cd(5);
    hSize_f2Hits_preL0->DrawClone("colz");
    cEnergy->cd(6);
    hSize_f2Hits_postL0->DrawClone("colz");
}