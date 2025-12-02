#include "ActKinematics.h"

#include "../PrettyStyle.C"

void DrawKinTheo()
{
    PrettyStyle(false);
    auto kin = ActPhysics::Kinematics("11Li(d,t)@82.5");
    kin.Draw();
}