#include "Rtypes.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TVirtualPad.h"
void pres()
{
    auto* gs {new TGraphErrors {"./Outputs/xs/g0_xs.dat", "%lg %lg %lg"}};
    gs->SetTitle("^{7}Li(d,d) g.s;#theta_{CM} [#circ];d#sigma/d#Omega [mb/sr]");
    gs->SetMarkerStyle(21);

    // Theoretical
    auto* gtheo {new TGraphErrors{"./Inputs/gsDA1p/fort.201", "%lg %lg"}};
    gtheo->SetLineWidth(2);
    gtheo->SetLineColor(kBlue);
    gtheo->Scale(1);

    auto* c0 {new TCanvas {"c0", "Presentation canvas"}};
    gPad->SetLogy();
    gs->SetMinimum(10);
    gs->SetMaximum(2e3);
    gs->Draw("ap");
    gtheo->Draw("c");
}
