#include "TROOT.h"
#include "TStyle.h"

void PrettyStyle(bool showStats = true)
{
    TStyle* st = new TStyle("PrettyStyle", "Pretty ROOT plots");

    // --- Colores y fuentes ---
    st->SetPalette(kViridis);
    st->SetTextFont(42);
    st->SetLabelFont(42, "XYZ");
    st->SetTitleFont(42, "XYZ");

    // --- Márgenes ---
    st->SetPadLeftMargin(0.13);
    st->SetPadRightMargin(0.10);
    st->SetPadBottomMargin(0.12);
    st->SetPadTopMargin(0.08);

    // --- Tamaños y offsets ---
    st->SetTitleSize(0.045, "XYZ");
    st->SetLabelSize(0.040, "XYZ");
    st->SetTitleOffset(1.10, "X");
    st->SetTitleOffset(1.25, "Y");

    // --- Cuadrícula ---
    st->SetPadGridX(true);
    st->SetPadGridY(true);

    // --- Líneas y marcadores ---
    st->SetHistLineWidth(2);
    st->SetFrameLineWidth(2);
    st->SetLineWidth(2);
    st->SetMarkerSize(1.2);

    // --- Fondo ---
    st->SetCanvasColor(0);
    st->SetPadColor(0);
    st->SetFrameFillColor(0);

    // --- Leyendas ---
    st->SetLegendBorderSize(0);
    st->SetLegendFillColor(0);

    // --- Estadísticas (activables/desactivables) ---
    if(showStats)
    {
        st->SetOptStat(1110); // entries, mean, rms
        st->SetOptFit(111);   // parámetros del fit
    }
    else
    {
        st->SetOptStat(0);
        st->SetOptFit(0);
    }

    // --- Título centrado (versión corregida) ---
    st->SetTitleAlign(23); // centered
    st->SetTitleX(0.5);    // center in X
    st->SetTitleY(0.93);   // slightly below top
    st->SetTitleW(0.30);   // width: 30% of canvas (no más!)
    st->SetTitleH(0.05);   // height

    // --- Aplicar estilo ---
    st->cd();
    gROOT->ForceStyle();
}
