#ifndef PRETTYSTYLE_H
#define PRETTYSTYLE_H

#include "TColor.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"

void PrettyStyle(bool showStats = true, bool palete2D = true)
{
    TStyle* st = new TStyle("PrettyStyle", "Pretty ROOT plots");

    // --- PALETA DE COLORES DISTINTIVA PARA LÍNEAS ---
    // Tableau 10 - colores bien diferenciados
    const Int_t NColors = 10;
    Int_t palette[NColors] = {
        TColor::GetColor("#1B3A6F"), // 0 - Azul oscuro   -> puntos
        TColor::GetColor("#6A3D9A"), // 5 - Violeta
        TColor::GetColor("#E66100"), // 6 - Naranja        (empieza a acercarse al rojo)
        TColor::GetColor("#C2185B"), // 7 - Magenta
        TColor::GetColor("#8B1E1E"), // 8 - Rojo oscuro    (tono casi igual al del fit)
        TColor::GetColor("#008C95"), // 1 - Cian oscuro    (muy distinto del rojo)
        TColor::GetColor("#008C45"), // 2 - Verde fuerte
        TColor::GetColor("#0057B8"), // 3 - Azul fuerte
        TColor::GetColor("#176B3A"), // 4 - Verde oscuro
        TColor::GetColor("#D62728"), // 9 - Rojo fuerte    -> fit total
    };
    // --- Colores y fuentes ---
    if(palete2D)
        st->SetPalette(kCividis);
    else
        st->SetPalette(NColors, palette);

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
    st->SetGridColor(17); // gris muy suave
    st->SetGridStyle(3);  // punteado fino
    st->SetGridWidth(1);

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
    st->SetLegendBorderSize(0);
    st->SetLegendFillColor(0);
    st->SetLegendFont(42);
    st->SetLegendTextSize(0.035);

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
        st->SetOptTitle(0);
    }

    // --- Título centrado ---
    st->SetTitleAlign(23);
    st->SetTitleX(0.5);
    st->SetTitleY(0.93);
    st->SetTitleW(0.30);
    st->SetTitleH(0.05);

    // st->SetOptTitle(0); // <--- apaga el título por completo

    // --- Ticks y divisiones ---
    st->SetPadTickX(1);
    st->SetPadTickY(1);
    st->SetNdivisions(510, "X");
    st->SetNdivisions(510, "Y");

    // --- Aplicar estilo ---
    st->cd();
    gROOT->ForceStyle();
}

#endif