#include "ActKinematics.h"

#include "TGraph.h"
#include "TAxis.h"
#include "TColor.h"
#include "../PrettyStyle.C"

void DrawKinTheo()
{
    PrettyStyle(false);

    // Kinematics
    auto kin11dd = ActPhysics::Kinematics("11Li(d,d)@82.5");
    auto kin11dp = ActPhysics::Kinematics("11Li(d,p)@82.5");
    auto kin7dd  = ActPhysics::Kinematics("7Li(d,d)@52.5");
    auto kin7dp  = ActPhysics::Kinematics("7Li(d,p)@52.5");

    auto* theo11dd = kin11dd.GetKinematicLine3();
    auto* theo11dp = kin11dp.GetKinematicLine3();
    auto* theo7dd  = kin7dd.GetKinematicLine3();
    auto* theo7dp  = kin7dp.GetKinematicLine3();

    // Colores
    int col_dd = kBlue+1;   // (d,d)
    int col_dp = kRed+1;    // (d,p)

    // Estilos de línea:
    theo11dd->SetLineColor(col_dd);
    theo11dd->SetLineStyle(1);
    theo11dd->SetLineWidth(3);

    theo11dp->SetLineColor(col_dp);
    theo11dp->SetLineStyle(1);
    theo11dp->SetLineWidth(3);

    theo7dd->SetLineColor(col_dd);
    theo7dd->SetLineStyle(10);
    theo7dd->SetLineWidth(3);

    theo7dp->SetLineColor(col_dp);
    theo7dp->SetLineStyle(10);
    theo7dp->SetLineWidth(3);

    // Dibujar
    theo11dd->Draw();
    theo11dd->GetXaxis()->SetLimits(0, 180);  // <-- límite X
    theo11dp->Draw("same");
    theo7dd->Draw("same");
    theo7dp->Draw("same");

    // Leyenda
    auto* leg = new TLegend(0.15, 0.70, 0.42, 0.92);  // posición (NDC)
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.055);

    leg->SetEntrySeparation(0.05);

    leg->AddEntry(theo11dd, "^{11}Li(d,d)", "l");
    leg->AddEntry(theo11dp, "^{11}Li(d,p)", "l");
    leg->AddEntry(theo7dd,  "^{7}Li(d,d)", "l");
    leg->AddEntry(theo7dp,  "^{7}Li(d,p)", "l");

    leg->Draw();
}
