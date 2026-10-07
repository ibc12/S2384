#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActParticle.h"
#include "ActRunner.h"
#include "ActSRIM.h"
#include "ActSilData.h"
#include "ActSilMatrix.h"
#include "ActSilSpecs.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TH2.h"
#include "TH2D.h"
#include "TMath.h"
#include "TRandom.h"
#include "TString.h"

#include "Math/AxisAngle.h"
#include "Math/DisplacementVector3D.h"
#include "Math/Point3Dfwd.h"
#include "Math/Rotation3D.h"
#include "Math/RotationZYX.h"
#include "Math/Vector3D.h"

#include <fstream>
#include <map>
#include <string>

#include "../../PostAnalysis/HistConfig.h"

using XYZPoint = ROOT::Math::XYZPointF;
using XYZVector = ROOT::Math::XYZVectorF;

double GetPhi3D(const XYZVector& beam, const XYZVector& other)
{
    auto beamFrame {beam.Unit()};
    XYZVector worldFrame {1, 0, 0};

    // Rotación inversa a la de la simulación: lleva el haz -> X (lab -> frame del haz)
    // Eje invertido (u x X en vez de X x u), mismo ángulo
    auto cross {beamFrame.Cross(worldFrame)};
    auto angle {TMath::ACos(beamFrame.Dot(worldFrame))};
    ROOT::Math::AxisAngle axis {cross, angle};
    ROOT::Math::Rotation3D rotation {axis};

    // Traza en el frame del haz
    auto t {rotation(other.Unit())};
    return TMath::ATan2(t.Y(), t.Z()) * TMath::RadToDeg();
}

double GetPhi3DLegacy(const XYZVector& beam, const XYZVector& other)
{
    // TODO: Check validity of phi calculation

    // auto ub {beam.Unit()};            // unitary beam
    auto trackUnitary {other.Unit()};
    // XYZVector yz {0, ub.Y(), ub.Z()}; // beam dir in YZ plane
    // auto dot {other.Unit().Dot(yz) / yz.R()};
    // return TMath::ACos(dot) * TMath::RadToDeg();
    return TMath::ATan2(trackUnitary.Y(), trackUnitary.Z()) * TMath::RadToDeg();
}

// Igual que ActSim::Runner::RotateToWorldFrame: frame del haz -> laboratorio
XYZVector RotateToWorldFrame(const XYZVector& vBeamFrame, const XYZVector& beamDir)
{
    const auto u {beamDir.Unit()};
    const XYZVector x {1, 0, 0};
    const auto axis {x.Cross(u)};
    const double s {axis.R()};
    ROOT::Math::Rotation3D rot; // identidad por defecto
    if(s >= 1e-12)
        rot = ROOT::Math::Rotation3D {ROOT::Math::AxisAngle {axis.Unit(), TMath::ATan2(s, x.Dot(u))}};
    return rot(vBeamFrame);
}

void EnsureGoodDefinitionPhi()
{
    // X = cos(th), Y = sin(th) sin(phi), Z = sin(th) cos(phi)
    auto dir = [](double thetaDeg, double phiDeg) -> XYZVector
    {
        const double th {thetaDeg * TMath::DegToRad()};
        const double ph {phiDeg * TMath::DegToRad()};
        return {static_cast<float>(TMath::Cos(th)), static_cast<float>(TMath::Sin(th) * TMath::Sin(ph)),
                static_cast<float>(TMath::Sin(th) * TMath::Cos(ph))};
    };

    const double phiBeam {5};   // dirección del haz en el laboratorio
    const double phiLight {45}; // phi de la ligera respecto al haz
    for(double thetaBeam : {5., 15., 45.})
    {
        const auto beam {dir(thetaBeam, phiBeam)};
        std::cout << "== beam: theta = " << thetaBeam << ", phi = " << phiBeam << " ==\n";
        for(double thetaOther : {20., 40., 60., 80.})
        {
            // Ligera generada en el frame del haz, y luego llevada al laboratorio
            const auto other {RotateToWorldFrame(dir(thetaOther, phiLight), beam)};
            std::cout << "  theta(vs haz) = " << thetaOther
                      << "  check = " << TMath::ACos(beam.Unit().Dot(other.Unit())) * TMath::RadToDeg()
                      << "  phi(new) = " << GetPhi3D(beam, other) << "  phi(legacy) = " << GetPhi3DLegacy(beam, other)
                      << '\n';
        }
    }
}