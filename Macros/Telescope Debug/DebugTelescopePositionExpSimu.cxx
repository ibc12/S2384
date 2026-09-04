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


void DebugTelescopePositionExpSimu()
{
    std::string beam {"11Li"};
    // ROOT::EnableImplicitMT();
    ROOT::RDataFrame dfini {"PreProcessed_Tree",
                            TString::Format("../../PostAnalysis/Outputs/tree_preprocess_F_%s.root", beam.c_str())};
    auto df {dfini.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == true; }, {"MergerData"})};

    // Create histograms for post L0 change, and for pre.
    // This are the post
    auto hSPf2 {HistConfig::SP.GetHistogram()};
    hSPf2->SetTitle("SP for f2");
    auto hSPf3 {HistConfig::SP.GetHistogram()};
    hSPf3->SetTitle("SP for f3");
    auto hAllSPf3 {HistConfig::SP.GetHistogram()};
    hAllSPf3->SetTitle("All SP projection to f3");
    auto hVertexOnf3 {HistConfig::SP.GetHistogram()};
    hVertexOnf3->SetTitle("Vertex projection on f3");
    // Now the pre
    auto hSPf2_pre {HistConfig::SP.GetHistogram()};
    hSPf2_pre->SetTitle("SP for f2 pre L0 change");
    auto hSPf3_pre {HistConfig::SP.GetHistogram()};
    hSPf3_pre->SetTitle("SP for f3 pre L0 change");
    auto hAllSPf3_pre {HistConfig::SP.GetHistogram()};
    hAllSPf3_pre->SetTitle("All SP projection to f3 pre L0 change");
    auto hVertexOnf3_pre {HistConfig::SP.GetHistogram()};
    hVertexOnf3_pre->SetTitle("Vertex projection on f3 pre L0 change");

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

    // For each event that hit f2 and f3 compute the silicon point and fill the histograms
    dfF2F3.Foreach(
        [&hSPf2, &hSPf3, &hVertexOnf3, &hSPf2_pre, &hSPf3_pre, &hVertexOnf3_pre, &runner, &sils,
         &sils_pre](ActRoot::SilData& sil, ActRoot::MergerData& m, ActRoot::TPCData& tpc)
        {
            int beamIdx = m.fBeamIdx;
            auto beamDirFloat = tpc.fClusters[beamIdx].GetLine().GetDirection().Unit();
            ROOT::Math::XYZVector beamDir {beamDirFloat.X(), beamDirFloat.Y(), beamDirFloat.Z()};

            int heavyIdx = m.fHeavyIdx;
            auto heavyDirFloat = tpc.fClusters[heavyIdx].GetLine().GetDirection().Unit();
            ROOT::Math::XYZVector heavyDir {heavyDirFloat.X(), heavyDirFloat.Y(), heavyDirFloat.Z()};

            ROOT::Math::XYZPoint vertex {m.fRP.X(), m.fRP.Y(), m.fRP.Z()};

            int silIndexHeavy {};
            int silIndexHeavyF2 {};
            ROOT::Math::XYZPoint silPointHeavy {};
            ROOT::Math::XYZPoint silPointHeavyF2 {};
            if(m.fRun > 66)
            {
                std::tie(silIndexHeavy, silPointHeavy) = sils->FindSPInLayer("f3", vertex, heavyDir);
                std::tie(silIndexHeavyF2, silPointHeavyF2) = sils->FindSPInLayer("f2", vertex, heavyDir);
                hSPf2->Fill(silPointHeavyF2.Y(), silPointHeavyF2.Z());
                hSPf3->Fill(silPointHeavy.Y(), silPointHeavy.Z());
                hVertexOnf3->Fill(vertex.Y(), vertex.Z());
            }
            else
            {
                std::tie(silIndexHeavy, silPointHeavy) = sils_pre->FindSPInLayer("f3", vertex, heavyDir);
                std::tie(silIndexHeavyF2, silPointHeavyF2) = sils_pre->FindSPInLayer("f2", vertex, heavyDir);
                hSPf2_pre->Fill(silPointHeavyF2.Y(), silPointHeavyF2.Z());
                hSPf3_pre->Fill(silPointHeavy.Y(), silPointHeavy.Z());
                hVertexOnf3_pre->Fill(vertex.Y(), vertex.Z());
            }
        },
        {"SilData", "MergerData", "TPCData"});

    // And for all events fill theextrapolated to f3
    df.Foreach(
        [&hAllSPf3, &hAllSPf3_pre, &sils, &sils_pre](ActRoot::SilData& sil, ActRoot::MergerData& m,
                                                     ActRoot::TPCData& tpc)
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
            if(m.fRun > 66)
            {
                extrapolatedPointF3 = heavyLine.MoveToX(sils->GetLayer("f3").GetPoint().X());
                hAllSPf3->Fill(extrapolatedPointF3.Y(), extrapolatedPointF3.Z());
            }
            else
            {
                extrapolatedPointF3 = heavyLine.MoveToX(sils_pre->GetLayer("f3").GetPoint().X());
                hAllSPf3_pre->Fill(extrapolatedPointF3.Y(), extrapolatedPointF3.Z());
            }
        },
        {"SilData", "MergerData", "TPCData"});

    // Save events in a .dat file
    // std::ofstream outfileF3("./Outputs/events_f3_only.dat");
    // dfF3.Foreach([&outfileF3](ActRoot::MergerData m) { m.Stream(outfileF3); }, {"MergerData"});

    // Plot
    auto* cSP {new TCanvas {"cSP", "Sil Points for post L0 change", 1200, 800}};
    cSP->Divide(2, 2);
    cSP->cd(1);
    hSPf2->DrawClone("colz");
    sils->GetLayer("f2").GetSilMatrix()->Draw();
    cSP->cd(2);
    hSPf3->DrawClone("colz");
    auto f3 = sils->GetLayer("f3").GetSilMatrix();
    f3->DrawClone();
    cSP->cd(3);
    hVertexOnf3->DrawClone("colz");
    f3->DrawClone();
    cSP->cd(4);
    hAllSPf3->DrawClone("colz");
    f3->DrawClone();

    auto* cSP_pre {new TCanvas {"cSP_pre", "Sil Points for pre L0 change", 1200, 800}};
    cSP_pre->Divide(2, 2);
    cSP_pre->cd(1);
    hSPf2_pre->DrawClone("colz");
    sils_pre->GetLayer("f2").GetSilMatrix()->Draw();
    cSP_pre->cd(2);
    hSPf3_pre->DrawClone("colz");
    auto f3_pre = sils_pre->GetLayer("f3").GetSilMatrix();
    f3_pre->DrawClone();
    cSP_pre->cd(3);
    hVertexOnf3_pre->DrawClone("colz");
    f3_pre->DrawClone();
    cSP_pre->cd(4);
    hAllSPf3_pre->DrawClone("colz");
    f3_pre->DrawClone();
}