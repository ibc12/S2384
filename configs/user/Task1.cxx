#include "Task1.h"

#include "ActMergerData.h"
#include "ActSRIM.h"
#include "ActSilData.h"
#include "ActVTask.h"

#include <iostream>

TSpline3* ActAlgorithm::Task1::BuildSRIMspline(ActPhysics::SRIM* srim, const std::string& particleKey, double range,
                                               double step, double sOffset)
{
    std::vector<double> s_pts, y_pts;

    double E = srim->EvalInverse(particleKey, range);
    bool ran_out = false;

    for(double r = 0; r < range; r += step)
    {
        double Epost = srim->Slow(particleKey, E, step);
        if(Epost < 0.0)
        {
            Epost = 0.0;
            ran_out = true;
        }

        double dE = E - Epost;
        E = Epost;

        if(dE <= 0)
        {
            if(ran_out)
                break;
            else
                continue;
        }

        s_pts.push_back(r + 0.5 * step + sOffset);
        y_pts.push_back(dE);

        if(ran_out)
            break;
    }

    // Ancla final en cero
    if(!y_pts.empty() && y_pts.back() > 0.0)
    {
        s_pts.push_back(range);
        y_pts.push_back(0.0);
    }

    int nSteps = static_cast<int>(s_pts.size());
    if(nSteps < 3)
    {
        std::cout << "[Task1] Not enough points for spline of '" << particleKey << "' (nSteps=" << nSteps << ")\n";
        return nullptr;
    }

    auto* sp = new TSpline3(("spSRIM_" + particleKey).c_str(), s_pts.data(), y_pts.data(), nSteps, "b1,e1", 0, 0);
    sp->SetNpx(3000);
    std::cout << "[Task1] Spline '" << particleKey << "' — xmax=" << sp->GetXmax()
              << "  val(range-1)=" << sp->Eval(range - 1) << "\n";
    return sp;
}

double FindPositionFromChargeFraction(TH1* h, double frac)
{
    double total = h->Integral("width");
    double accum = 0.0;

    for(int i = 1; i <= h->GetNbinsX(); ++i)
    {
        accum += h->GetBinContent(i) * h->GetBinWidth(i);

        if(accum / total >= frac)
            return h->GetBinCenter(i);
    }

    return h->GetXaxis()->GetXmax();
}

TF1* FitSRIMtoChargeProfileFixedEnd(TH1* hCharge, TSpline3* spSRIM, const std::string& particleKey,
                                    double maxEndShift = 20.0) // cuanto puede moverse el final (mm)
{
    if(!spSRIM)
        return nullptr;

    // --- región útil del perfil ---
    double sOffset = FindPositionFromChargeFraction(hCharge, 0.02);
    double sEndData = FindPositionFromChargeFraction(hCharge, 0.98);
    double sEndFull = FindPositionFromChargeFraction(hCharge, 0.9999);

    std::cout << "Fit start offset = " << sOffset << " mm\n";
    std::cout << "Nominal end position = " << sEndData << " mm\n";

    // --- final físico del spline (Bragg peak final) ---
    double rMax = spSRIM->GetXmax();

    std::string fname = "fSRIMfitFixedEnd_" + particleKey;

    TF1* f = new TF1(
        fname.c_str(),
        [spSRIM, sOffset, sEndData, rMax](double* x, double* par)
        {
            double A = par[0];
            double deltaEnd = par[1];

            double s = x[0];

            if(s < sOffset)
                return 1e-9;

            // posición efectiva del final del track
            double sEndFit = sEndData + deltaEnd;

            // alineación spline-datos
            double r = rMax - (sEndFit - s);

            if(r < spSRIM->GetXmin() || r > spSRIM->GetXmax())
                return 1e-9;

            return A * spSRIM->Eval(r);
        },
        sOffset, sEndData + maxEndShift, 2);

    // --- inicialización ---
    double integral = hCharge->Integral("width");
    double initAmp = integral / hCharge->GetNbinsX();

    f->SetParameters(initAmp, 0.0);

    f->SetParName(0, "Amplitude");
    f->SetParName(1, "EndShift");

    f->SetParLimits(0, 0, 1e12);
    f->SetParLimits(1, -maxEndShift, maxEndShift);

    f->SetNpx(3000);

    hCharge->Fit(f, "QR0", "", sOffset, sEndData + maxEndShift);

    return f;
}

bool ActAlgorithm::Task1::Run()
{
    // Get the charge profile and fit it to get a range value
    if(!fSplinesBuilt)
    {
        auto* srim = new ActPhysics::SRIM;
        srim->ReadTable("d", "../../Calibrations/SRIM/2H_900mb_CF4_95-5.txt");
        fSpline = BuildSRIMspline(srim, "d", 250, 0.5);

        fSplinesBuilt = true;
    }


    return true;
}

void ActAlgorithm::Task1::Print()
{
    std::cout << "Task1 prints!" << '\n';
}

// Create symbol to load class from .so
// extern "C" disables C++ function name mangling
extern "C" ActAlgorithm::Task1* Create()
{
    return new ActAlgorithm::Task1;
}