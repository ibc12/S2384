#include "ActDataManager.h"
#include "ActSilData.h"
#include "ActTPCData.h"
#include "ActModularData.h"
#include "ActModularData.h"
#include "ActMergerData.h"
#include "ActTPCData.h"
#include "ActCluster.h"
#include "ActVoxel.h"
#include "ActCutsManager.h"

#include <ROOT/RDataFrame.hxx>

#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"
#include "TF1.h"
#include "TLatex.h"
#include "TLine.h"
#include "TEllipse.h"

void vDriftfromAlpha()
{
    ROOT::EnableImplicitMT();
    // Get the data run 123
    auto df{ROOT::RDataFrame("GETTree", "../../RootFiles/Cluster/Clusters_Run_0123.root")};

    // df.Describe().Print();

    // Define last point of cluster in x y z, as the projection of the alpha track
    auto dfLastPoint = df.Define("fLastPoint", [](ActRoot::TPCData &d)
                                 {
        if(d.fClusters.size() != 1)
        {
            ROOT::Math::XYZPointF lastPoint {-1000, -1000, -1000};
            return lastPoint;
        }
        else
        {
            auto cluster {d.fClusters[0]};
            auto line {cluster.GetRefToLine()};
            auto dir {line.GetDirection()};
            cluster.SortAlongDir(dir);
            auto lastVoxel {cluster.GetRefToVoxels().back()};
            auto projectionPointLine {line.ProjectionPointOnLine(lastVoxel.GetPosition())};
            return projectionPointLine;
        } }, {"TPCData"});

    // Define other point
    auto dfOtherPoint = dfLastPoint.Define("fOtherPoint", [](ActRoot::TPCData &d)
                                           {
        if(d.fClusters.size() != 1)
        {
            ROOT::Math::XYZPointF otherPoint {-1000, -1000, -2000};
            return otherPoint;
        }
        else
        {
            auto cluster {d.fClusters[0]};
            auto line {cluster.GetRefToLine()};
            auto otherPoint {line.MoveToX(-50)};
            return otherPoint;
        } }, {"TPCData"});

    // Define the coordinates of the points
    auto dfXY = dfOtherPoint
                    .Define("fLastX", "fLastPoint.X()")
                    .Define("fLastY", "fLastPoint.Y()")
                    .Define("fOtherX", "fOtherPoint.X()")
                    .Define("fOtherY", "fOtherPoint.Y()");

    // Plot
    auto c = new TCanvas("c", "Points XY", 1000, 500);
    auto hLast = dfXY.Histo2D({"hLast", "LastPoint XY;X [pads];Y [pads]", 1000, -200, 200, 1000, -200, 200}, "fLastX", "fLastY");
    auto hOther = dfXY.Histo2D({"hOther", "OtherPoint XY;X [pads];Y [pads]", 1000, -200, 200, 1000, -200, 200}, "fOtherX", "fOtherY");
    hLast->DrawClone("colz");
    hOther->DrawClone("same");

    // Create the lines and draw them with the foreach
    int counter = 0;
    dfXY.Foreach(
        [&](float otherX, float otherY, float lastX, float lastY)
        {
            counter++;
            auto line = new TLine(otherX, otherY, lastX, lastY);
            line->SetLineColorAlpha(kBlue, 0.3); // transparente para ver cruces
            if (counter % 50 == 0 && lastX > 5 && lastX < 60)
                line->Draw("same");
        },
        {"fOtherX", "fOtherY", "fLastX", "fLastY"});

    // With the plot I guess the alpha source is at (-27, 41)
    float xSource{-27.};
    float ySource{41.};

    auto dfDrift = dfXY.Define("fDeltaZ", [&](ActRoot::TPCData &d)
                               {
        if(d.fClusters.size() != 1)
            return -1000.;
        else
        {
            auto cluster {d.fClusters[0]};
            auto line {cluster.GetRefToLine()};
            auto dir {line.GetDirection()};
            cluster.SortAlongDir(dir);
            //auto firstVoxel {cluster.GetRefToVoxels().front()};
            //auto projectionFirstPointLine {line.ProjectionPointOnLine(firstVoxel.GetPosition())};
            auto lastVoxel {cluster.GetRefToVoxels().back()};
            auto projectionLastPointLine {line.ProjectionPointOnLine(lastVoxel.GetPosition())};
            auto zSource {line.MoveToX(xSource).Z()};
            double deltaZ = projectionLastPointLine.Z() - zSource;
            return deltaZ *0.32; // Conversion factor from btb to micro seconds
        } }, {"TPCData"})
                       .Define("fLxy", [&](ActRoot::TPCData &d)
                               {
        if(d.fClusters.size() != 1)
            return -1000.;
        else
        {
            auto cluster {d.fClusters[0]};
            auto line {cluster.GetRefToLine()};
            auto dir {line.GetDirection()};
            cluster.SortAlongDir(dir);
            auto lastVoxel {cluster.GetRefToVoxels().back()};
            auto projectionPointLine {line.ProjectionPointOnLine(lastVoxel.GetPosition())};
            double lxy = TMath::Sqrt(TMath::Power(projectionPointLine.X() - xSource, 2) + TMath::Power(projectionPointLine.Y() - ySource, 2));
            return lxy * 2; // Conversion factor from pads to mm
        } }, {"TPCData"})
                       .Define("fDeltaZSquare", "fDeltaZ * fDeltaZ")
                       .Define("fLxySquare", "fLxy * fLxy");

    // Plot the DeltaZ and lxy
    auto graphDrift = dfDrift.Graph("fDeltaZ", "fLxy");
    graphDrift->SetTitle("Delta Z vs Lxy;#Delta Z [#mus]; Lxy [mm]");

    // Linearize the graph
    auto graphDriftLinear = dfDrift.Graph("fDeltaZSquare", "fLxySquare");
    graphDriftLinear->SetTitle("Delta Z^2 vs Lxy^2;#Delta Z^2 [#mus^2]; Lxy^2 [mm^2]");


    auto c1 = new TCanvas("c1", "Delta Z vs Lxy", 1400, 800);
    c1->DivideSquare(2);
    c1->cd(1);
    graphDrift->DrawClone("AP");
    c1->cd(2);
    graphDriftLinear->DrawClone("AP");

    // Cuts for good events (no broad region) and for each line
    ActRoot::CutsManager<std::string> cuts;
    // Gas PID
    cuts.ReadCut("goodEvents", "./Inputs/cut_DriftVelocity_GoodAlphaEvents.root");
    cuts.ReadCut("first", "./Inputs/cut_firstPeak.root");
    cuts.ReadCut("second", "./Inputs/cut_secondPeak.root");
    cuts.ReadCut("third", "./Inputs/cut_thirdPeak.root");

    auto dfFiltered = dfDrift.Filter([&](double lxy, double deltaZ)
                                { return cuts.IsInside("goodEvents", lxy, deltaZ); },
                                {"fLxy", "fDeltaZ"});

    auto graphDriftFiltered = dfFiltered.Graph("fDeltaZ", "fLxy");
    graphDriftFiltered->SetTitle("Delta Z vs Lxy ;#Delta Z [#mus]; Lxy [mm]");
    auto graphDriftFilteredLinear = dfFiltered.Graph("fDeltaZSquare", "fLxySquare");
    graphDriftFilteredLinear->SetTitle("Delta Z^2 vs Lxy^2;#Delta Z^2 [#mus^2];Lxy^2 [mm^2]");
    auto c2 = new TCanvas("c2", "Delta Z vs Lxy filtered", 1400, 800);
    c2->DivideSquare(2);
    c2->cd(1);
    graphDriftFiltered->DrawClone("AP");
    c2->cd(2);
    graphDriftFilteredLinear->DrawClone("AP");

    // Do graphs for each peak
    auto dfFirst = dfFiltered.Filter([&](double lxy2, double deltaZ2)
                                { return cuts.IsInside("first", lxy2, deltaZ2); },
                                {"fLxySquare", "fDeltaZSquare"});
    auto graphDriftLineFirst = dfFirst.Graph("fDeltaZSquare", "fLxySquare");
    graphDriftLineFirst->SetTitle("Delta Z^2 vs Lxy^2 (first peak);#Delta Z^2 [#mus^2];Lxy^2 [mm^2]");
    graphDriftLineFirst->Fit("pol1");
    auto f1 {graphDriftLineFirst->GetFunction("pol1")};
    f1->SetLineColor(kRed);
    auto dfSecond = dfFiltered.Filter([&](double lxy2, double deltaZ2)
                                { return cuts.IsInside("second", lxy2, deltaZ2); },
                                {"fLxySquare", "fDeltaZSquare"});
    auto graphDriftLineSecond = dfSecond.Graph("fDeltaZSquare", "fLxySquare");
    graphDriftLineSecond->SetTitle("Delta Z^2 vs Lxy^2 (second peak);#Delta Z^2 [#mus^2];Lxy^2 [mm^2]");
    graphDriftLineSecond->Fit("pol1");
    auto f2 {graphDriftLineSecond->GetFunction("pol1")};
    auto dfThird = dfFiltered.Filter([&](double lxy2, double deltaZ2)
                                { return cuts.IsInside("third", lxy2, deltaZ2); },
                                {"fLxySquare", "fDeltaZSquare"});
    auto graphDriftLineThird = dfThird.Graph("fDeltaZSquare", "fLxySquare");
    graphDriftLineThird->SetTitle("Delta Z^2 vs Lxy^2 (third peak);#Delta Z^2 [#mus^2];Lxy^2 [mm^2]");
    graphDriftLineThird->Fit("pol1");
    auto f3 {graphDriftLineThird->GetFunction("pol1")};
    auto c3 = new TCanvas("c3", "Delta Z vs Lxy lines", 2100, 700);
    c3->DivideSquare(3);
    c3->cd(1);
    graphDriftLineFirst->DrawClone("AP");
    f1->Draw("same");
    c3->cd(2);
    graphDriftLineSecond->DrawClone("AP");
    f2->Draw("same");
    c3->cd(3);
    graphDriftLineThird->DrawClone("AP");
    f3->Draw("same");

    // Draw them also in the filtered plot
    c2->cd(2);
    f1->DrawClone("same");
    f2->SetLineColor(kGreen);
    f2->DrawClone("same");
    f3->SetLineColor(kBlue);
    f3->DrawClone("same");
    // Text of the fit parameters
    auto t1 = new TLatex(100, 23000, TString::Format("First peak: Vdrift = %.2f #pm ", TMath::Sqrt(-f1->GetParameter(1))));
    auto t2 = new TLatex(100, 21000, TString::Format("Second peak: Vdrift = %.2f #pm ", TMath::Sqrt(-f2->GetParameter(1))));
    auto t3 = new TLatex(100, 19000, TString::Format("Third peak: Vdrift = %.2f #pm ", TMath::Sqrt(-f3->GetParameter(1))));
    t1->DrawClone();
    t2->DrawClone();
    t3->DrawClone();
    // Add in the filtered plot the elipses
    // TEllipse(x_center, y_center, r_x, r_y, theta1, theta2, angle_rotation)
    // theta1/theta2 = ángulos de barrido (0–360), angle_rotation = rotación de la elipse en grados
    c2->cd(1);
    // First peak ellipses
    auto e1_1st = new TEllipse(0.0, 0.0, 13.0, 116.0, 60.0, 120.0, 0.0); // rotada 30°
    e1_1st->SetFillStyle(0);
    e1_1st->SetNoEdges();
    auto e2_1st = new TEllipse(0.0, 0.0, 14.7, 124.0, 60.0, 120.0, 0.0); // rotada 30°
    e2_1st->SetFillStyle(0);
    e2_1st->SetNoEdges();
    e1_1st->SetLineColor(kRed);
    e1_1st->DrawClone("same");
    e2_1st->SetLineColor(kRed);
    e2_1st->DrawClone("same");
    // Second peak ellipses
    auto e1_2nd = new TEllipse(0.0, 0.0, 16.0, 128.0, 60.0, 120.0, 0.0); // rotada 30°
    e1_2nd->SetFillStyle(0);
    e1_2nd->SetNoEdges();
    auto e2_2nd = new TEllipse(0.0, 0.0, 17, 137.0, 60.0, 120.0, 0.0); // rotada 30°
    e2_2nd->SetFillStyle(0);
    e2_2nd->SetNoEdges();
    e1_2nd->SetLineColor(kGreen);
    e1_2nd->DrawClone("same");
    e2_2nd->SetLineColor(kGreen);
    e2_2nd->DrawClone("same");
    // Third peak ellipses
    auto e1_3rd = new TEllipse(0.0, 0.0, 17.2, 143.0, 60.0, 120.0, 0.0); // rotada 30°
    e1_3rd->SetFillStyle(0);
    e1_3rd->SetNoEdges();
    auto e2_3rd = new TEllipse(0.0, 0.0, 18.5, 151.0, 60.0, 120.0, 0.0); // rotada 30°
    e2_3rd->SetFillStyle(0);
    e2_3rd->SetNoEdges();
    e1_3rd->SetLineColor(kBlue);
    e1_3rd->DrawClone("same");
    e2_3rd->SetLineColor(kBlue);
    e2_3rd->DrawClone("same");

    
}
