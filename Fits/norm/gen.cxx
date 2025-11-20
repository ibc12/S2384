#include "PhysExperiment.h"

// Generates the different normalizations
void gen()
{
    // Number of incoming beams
    double Ntrigger_7Li {648000};
    double Ntrigger_11Li {3125130};
    double Ndiv {300};
    // Total ACTAR length used in the LISE calculation
    double totalLength {25.6}; // cm
    double Nd {4.126e19 * totalLength}; // number of deuterons in the target

    // Build
    PhysUtils::Experiment norm11Li {Nd, Ntrigger_11Li, Ndiv};
    norm11Li.Print();
    norm11Li.Write("./11Li_norm.dat");

    PhysUtils::Experiment norm7Li {Nd, Ntrigger_7Li, Ndiv};
    norm7Li.Print();
    norm7Li.Write("./7Li_norm.dat");
}