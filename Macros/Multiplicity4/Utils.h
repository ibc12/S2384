#ifndef UTILS_H
#define UTILS_H

#include "Math/Point3Dfwd.h"
#include "Math/Vector3D.h"


namespace Utils {

constexpr float scaleXY = 2.0f;
constexpr float scaleZ  = 2.84032f;

inline void ScalePoint(ROOT::Math::XYZPointF& point, bool addOffset = false)
{
    if(addOffset)
        point += ROOT::Math::XYZVector{0.5f, 0.5f, 0.5f};

    point.SetX(point.X() * scaleXY);
    point.SetY(point.Y() * scaleXY);
    point.SetZ(point.Z() * scaleZ);
}

double GetTheta3D(const ROOT::Math::XYZVectorF& beam, const ROOT::Math::XYZVectorF& other)
{
    auto dot {beam.Unit().Dot(other.Unit())};
    return TMath::ACos(dot) * TMath::RadToDeg();
}
} // namespace Utils

#endif