#ifndef Histos_h
#define Histos_h
#include "ROOT/RDF/HistoModels.hxx"
namespace Histos
{
// Ex
const ROOT::RDF::TH1DModel Ex {"hEx", "Ex;Ex [MeV];Counts", 400, -5, 10};

// Ecm
const ROOT::RDF::TH1DModel Ecm {"hEcm", "Ecm;E_{CM} [MeV];Counts", 300, -5, 10};

// T4
const ROOT::RDF::TH1DModel T4Lab {"hT4Lab", "Kinetic Energy of heavy;T4_{Lab} [MeV];Counts", 300, 0, 100};

// T3
const ROOT::RDF::TH1DModel T3Lab {"hT3Lab", "Kinetic Energy of light;T3_{Lab} [MeV];Counts", 100, 0, 20};

// T1
const ROOT::RDF::TH1DModel T1Lab {"hT1Lab", "Kinetic Energy of beam;T1_{Lab} [MeV];Counts", 100, 80, 85};

// Range In Gas 
const ROOT::RDF::TH1DModel RangeInGas {"hRangeInGas", "Range;Range [mm];Counts", 300, 0, 2000};

// Range In Gas * sin(theta)
const ROOT::RDF::TH2DModel RangeInGasSin {"hRangeInGasSin", "Range * sin(theta); #theta_{Lab} [deg];Rsin#theta", 180, 0, 180, 300, 0, 300};

// Range
const ROOT::RDF::TH1DModel RangeInSil {"hRangeInSil", "Range;Range [mm];Counts", 300, 0, 1};

// SP
const ROOT::RDF::TH2DModel SP {"hSP", "SP;X or Y [mm];Z [mm]", 200, -10, 300, 200, -10, 300};

// RP
const ROOT::RDF::TH2DModel RP {"hRP", "RP;X [mm];Y [mm]", 200, -10, 300, 200, -10, 300};

// RP ZY
const ROOT::RDF::TH2DModel RP_ZY {"hRP_ZY", "RP;Y [mm];Z [mm]", 200, -10, 300, 200, -10, 300};

// Distance to SP
const ROOT::RDF::TH1DModel DistanceSP {"hDistanceSP", "Distance to SP;Distance [mm]", 200, 0, 400};

// RP vs E
const ROOT::RDF::TH2DModel RP_E {"hRPE", "RPvsE;RP.X() [mm];E [MeV]", 200, -10, 300, 200, 0, 40};

// Eficiency RP
const ROOT::RDF::TH1DModel RP_eff {"hRPeff", "RPeff;RP.X() [mm]", 13, 0, 260};

// Straggling
const ROOT::RDF::TH1DModel Straggling {"hStraggling", "Straggling;Straggling [#mum]", 30, 0, 25};

// Esil vs EbeforeSil
const ROOT::RDF::TH2DModel EsilAftervsBefore {"hEsilAftervsBefore", "EsilAftervsBefore;E_{Sil} [MeV];#E_{before sil} [MeV]", 350, 0, 90, 350, 0,
                                70};

// PID heavy
const ROOT::RDF::TH2DModel PIDHeavy {"hPID", "PID;E_{Sil} [MeV];#DeltaE_{gas} [MeV]", 350, 0, 90, 350, 0,
                                30};

// PID heavy telescope
const ROOT::RDF::TH2DModel PIDHeavyTelescope {"hPIDTelescope", "PID;#DeltaE_{Sil0} [MeV];#DeltaE_{sil1} [MeV]", 350, 0, 30, 350, 0,
                                90};

// PID light
const ROOT::RDF::TH2DModel PIDLight {"hPID", "PID;E_{Sil} [MeV];#DeltaE_{gas} [MeV]", 350, 0, 50, 350, 0,
                                10};

// PID heavy length
const ROOT::RDF::TH2DModel PIDHeavylength {"hPIDlength", "PID length;E_{Sil} [MeV];#DeltaE_{gas} [MeV]", 350, 0, 90, 350, 0,
                                0.3};

// PID light length
const ROOT::RDF::TH2DModel PIDLightlength {"hPIDlength", "PID length;E_{Sil} [MeV];#DeltaE_{gas} [MeV]", 350, 0, 90, 350, 0,
                                0.03};

// Kin
const ROOT::RDF::TH2DModel Kin {"hKin", "Kinematics;#theta_{light, Lab} [#circ];E_{light} [MeV]", 350, 0, 165, 350, 0,
                                40};

// Kin heavy
const ROOT::RDF::TH2DModel KinHeavy {"hKinHeavy", "Kinematics;#theta_{heavy, Lab} [#circ];E_{heavy} [MeV]", 100, 0, 20, 350, 20,
                                100};
                                
// Theta Lab vs Theta CM 
const ROOT::RDF::TH2DModel ThetaCMThetaLab {"hThetaCMThetaLab", "ThetaCM vs ThetaLab; #theta_{CM} [deg]; #theta_{Lab} [deg]",
                                180, 0, 180, 180, 0, 180};

const ROOT::RDF::TH2DModel ThetaLabVertex {"hThetaLabVertex", "ThetaLab vs Vertex; #theta_{lab} [deg]; vertex.X() [mm]",
                                180, 0, 180, 256, 0, 256};

// Efficiency
const ROOT::RDF::TH1DModel ThetaCM {"hThetaCM", "ThetaCM;#theta_{CM} [#circ]", 600, 0, 180};
const ROOT::RDF::TH1DModel ThetaLab {"hThetaLab", "ThetaLab;#theta_{lab} [#circ]", 600, 0, 180};
const ROOT::RDF::TH1DModel ThetaLabHeavy {"hThetaLab", "ThetaLabHeavy;#theta_{lab} [#circ]", 200, 0, 10};
const ROOT::RDF::TH1DModel PhiLab {"hPhiLab", "PhiLab;#phi_{lab} [#circ]", 600, 0, 360};
} // namespace Histos

#endif // !Histos_h
