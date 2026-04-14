#ifndef Task1_h
#define Task1_h

#include "ActSRIM.h"
#include "ActVTask.h"

#include "TF1.h"
#include "TSpline.h"

#include <map>
#include <string>
#include <vector>

namespace ActAlgorithm
{
class Task1 : public VTask
{
public:
    Task1() : VTask("Task1") {}

    bool Run() override;
    void Print() override;

private:
    TSpline3* fSpline;
    bool fSplinesBuilt = false;

    TSpline3* BuildSRIMspline(ActPhysics::SRIM* srim, const std::string& particleKey = "d", double range = 250,
                              double step = 0.5, double sOffset = 0.0);
    TF1* FitSRIMtoChargeProfileFixedEnd(TH1* hCharge, TSpline3* spSRIM, const std::string& particleKey = "d",
                                        double sOffset = 0.0, double sEndData = 200.0, double maxEndShift = 20.0);
    double FindPositionFromChargeFraction(TH1* h, double frac);
};
} // namespace ActAlgorithm

#endif // !Task1_h