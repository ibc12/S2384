// ============================================================
// CheckLastPosYInsideExclusionZone.cxx
// Coge los eventos de SimulationTTree con LastPosY dentro de la zona de exclusion y calcula:
//   - si la LIGERA (3), con su rango y angulo, deberia salir de la zona (Y predicha al final de la traza)
//   - si la PESADA (4) sale de la zona y cuantos pads fuera podria aportar
//     (se supone que la pesada tiene energia de sobra para atravesar la camara)
// Uso:  root -l
//       .L CheckLastPosYInsideExclusionZone.cxx+
//       CheckLastPosYInsideExclusionZone()
// ============================================================
#include "ActSRIM.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TTree.h"

#include <cmath>
#include <iostream>
#include <set>
#include <string>

namespace ExZone
{
// ---- Geometria (igual que en la simulacion) ----
constexpr double pad = 2.0;   // mm
constexpr int nPadsTPC = 128; // pads por lado (256 mm)
constexpr int iyMin = 56;     // zona de exclusion en pads
constexpr int iyMax = 71;
constexpr double yLow = iyMin * pad;        // 112 mm
constexpr double yHigh = (iyMax + 1) * pad; // 144 mm
constexpr int nPadsThreshold = 19;
constexpr double stepLen = 0.5;    // mm, paso para recorrer la traza
constexpr double rangeHeavy = 1e4; // mm, la pesada atraviesa la camara (solo la limita la TPC)

// ---- Unidades de las ramas del arbol (cambiar si no coinciden) ----
constexpr bool theta4InDeg = true; // theta4Lab guardado en grados
constexpr bool phi4InDeg = false;  // phi4Lab guardado en radianes (como phi3CM)

inline bool InZone(double y)
{
    return y >= yLow && y < yHigh;
}

struct Path
{
    double xEnd {};
    double yEnd {};
    std::set<std::pair<int, int>> pads; // pads (ix, iy) FUERA de la zona tocados por la traza
};

// Recorre la traza en 3D (dx, dy = componentes X, Y del vector unitario) hasta len o hasta salir de la TPC
Path Propagate(double x0, double y0, double dx, double dy, double len)
{
    Path p;
    p.xEnd = x0;
    p.yEnd = y0;
    for(double s = 0; s <= len; s += stepLen)
    {
        double x = x0 + s * dx;
        double y = y0 + s * dy;
        int ix = int(std::floor(x / pad));
        int iy = int(std::floor(y / pad));
        if(ix < 0 || ix >= nPadsTPC || iy < 0 || iy >= nPadsTPC)
            break; // sale de la TPC
        p.xEnd = x;
        p.yEnd = y;
        if(iy < iyMin || iy > iyMax)
            p.pads.insert({ix, iy});
    }
    return p;
}
} // namespace ExZone

void CheckLastPosYInsideExclusionZone()
{
    std::string beam {"7Li"};
    std::string target {"2H"};
    std::string light {"2H"};
    std::string fileName {"../../Simulation/Outputs/" + beam + "/test_rpCut/" + target + "_" + light +
                          "_TRIUMF_Eex_0.000_nPS_0_pPS_0_L1_240rpCut_paralelization.root"};
    // phi lab de la ligera = phi3CM + phiShiftLight.
    //   0  -> caso sin XS (phi3CM = phi3Lab)
    //   pi -> caso con XS (en la simulacion phi3CM se guarda como phi - pi)
    // Si el histograma hDY sale desplazado / con el signo cambiado, prueba el otro valor.
    const double phiShiftLight {0.};
    using namespace ExZone;

    // ---- Entrada ----
    auto* fin = TFile::Open(fileName.c_str());
    if(!fin || fin->IsZombie())
    {
        std::cout << "No se puede abrir " << fileName << '\n';
        return;
    }
    auto* tree = fin->Get<TTree>("SimulationTTree");
    double lastPosY {}, vertexX {}, vertexY {}, EVertex {}, theta3Lab {}, phi3CM {}, theta4Lab {}, phi4Lab {},
        nPadsSim {}, TL {};
    tree->SetBranchAddress("LastPosY", &lastPosY);
    tree->SetBranchAddress("vertexX", &vertexX);
    tree->SetBranchAddress("vertexY", &vertexY);
    tree->SetBranchAddress("EVertex", &EVertex);
    tree->SetBranchAddress("theta3Lab", &theta3Lab); // grados
    tree->SetBranchAddress("phi3CM", &phi3CM);       // radianes
    tree->SetBranchAddress("theta4Lab", &theta4Lab);
    tree->SetBranchAddress("phi4Lab", &phi4Lab);
    tree->SetBranchAddress("nPads", &nPadsSim);
    tree->SetBranchAddress("TL", &TL);

    // ---- SRIM (rango de la ligera a partir de EVertex) ----
    ActPhysics::SRIM srim;
    const std::string path {"../Calibrations/SRIM/"};
    const std::string gas {"900mb_CF4_95-5"};
    srim.ReadTable("light", path + light + "_" + gas + ".txt");

    // ---- Histogramas ----
    auto* hDY = new TH1D("hDY", "Y_{end}^{ligera, pred} - LastPosY;#DeltaY [mm];Counts", 100, -50, 50);
    auto* hYend = new TH2D("hYend", "Ligera: Y final predicha vs LastPosY;LastPosY [mm];Y_{end} pred [mm]", 44, 100,
                           156, 100, 60, 200);
    auto* hPadsH = new TH1D("hPadsH", "Pads fuera de la zona (pesada, geometrico);Pads;Counts", 128, 0, 128);
    auto* hPadsU =
        new TH2D("hPadsU", "Pads: simulacion vs geometrico (L+H);nPads sim;nPads geo (L+H)", 60, 0, 120, 60, 0, 120);

    // ---- Contadores ----
    long nSel {}, nLightExits {}, nHeavyExits {}, nHeavyRescue {}, nVtxOut {};

    const Long64_t N = tree->GetEntries();
    for(Long64_t i = 0; i < N; i++)
    {
        tree->GetEntry(i);
        if(!InZone(lastPosY))
            continue;
        nSel++;
        if(!InZone(vertexY))
            nVtxOut++;

        // ---------- LIGERA ----------
        const double th3 = theta3Lab * TMath::DegToRad();
        const double ph3 = phi3CM + phiShiftLight;
        const double RL = (EVertex > 0) ? srim.EvalRange("light", EVertex) : TL;
        auto pL = Propagate(vertexX, vertexY, std::cos(th3), std::sin(th3) * std::sin(ph3), RL);

        // ---------- PESADA (atraviesa la camara) ----------
        const double th4 = theta4InDeg ? theta4Lab * TMath::DegToRad() : theta4Lab;
        const double ph4 = phi4InDeg ? phi4Lab * TMath::DegToRad() : phi4Lab;
        auto pH = Propagate(vertexX, vertexY, std::cos(th4), std::sin(th4) * std::sin(ph4), rangeHeavy);

        // ---------- Union de pads (como hace PadsOutExclusionZone) ----------
        std::set<std::pair<int, int>> padsU {pL.pads};
        padsU.insert(pH.pads.begin(), pH.pads.end());

        const bool lightExits = !pL.pads.empty();
        const bool heavyExits = !pH.pads.empty();
        if(lightExits)
            nLightExits++;
        if(heavyExits)
            nHeavyExits++;
        if(int(pL.pads.size()) < nPadsThreshold && int(padsU.size()) >= nPadsThreshold)
            nHeavyRescue++;

        hDY->Fill(pL.yEnd - lastPosY);
        hYend->Fill(lastPosY, pL.yEnd);
        hPadsH->Fill(pH.pads.size());
        hPadsU->Fill(nPadsSim, padsU.size());
    }

    std::cout << "\n==== Eventos con LastPosY en [" << yLow << ", " << yHigh << ") mm ====\n"
              << "Seleccionados:                              " << nSel << " de " << N << '\n'
              << "  vertexY fuera de la zona:                 " << nVtxOut << '\n'
              << "  Ligera sale de la zona (teorico):         " << nLightExits << '\n'
              << "  Pesada sale de la zona (teorico):         " << nHeavyExits << '\n'
              << "  Ligera < " << nPadsThreshold << " pads pero L+H >= " << nPadsThreshold
              << " (salvados por pesada): " << nHeavyRescue << '\n';

    auto* c = new TCanvas("cExcl", "Exclusion zone check");
    c->Divide(2, 2);
    c->cd(1);
    hDY->DrawClone();
    c->cd(2);
    hYend->DrawClone("colz");
    c->cd(3);
    hPadsH->DrawClone();
    c->cd(4);
    hPadsU->DrawClone("colz");

    fin->Close();
}