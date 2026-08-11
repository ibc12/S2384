#include "ActKinematics.h"
#include "ActParticle.h"

#include "PhysOMP.h"

#include <iostream>

void GetOMPsParameters()
{
    int Z = 3;           // Lithium
    int A_dd = 7;        // Lithium-7
    int A_dp = A_dd + 1; // Lithium-8
    double energy {};
    if(A_dd == 7)
    {
        energy = 7.26557; // 7Li beam time
    }
    else if(A_dd == 11)
    {
        energy = 7.37188; // 11Li beam time
    }

    // Get particles for the kinematics
    std::string beam = "7Li";
    std::string target = "d";
    std::string light = "p";
    std::string heavy = "8Li";
    ActPhysics::Particle pb {beam};
    ActPhysics::Particle pt {target};
    ActPhysics::Particle pl {light};
    ActPhysics::Particle ph {heavy};

    double Ex {2.255}; // Excitation energy of the heavy particle (8Li)

    auto kin {ActPhysics::Kinematics(pt, pb, ph, pl, energy * 2, Ex)};
    auto equivalentEnergy {kin.ComputeEquivalentBeamEnergy()}; // Enegy for the protons potential (as incident proton)

    std::cout << "Equivalent beam energy for Z = " << Z << ", A = " << A_dd << ", deuteron energy = " << energy * 2
              << " MeV is: " << equivalentEnergy << " MeV" << std::endl;


    int potIdx {0}; // index of the potential to use, if we want FRESCO to have multiple potentials (e.g., for different
                    // energy ranges)

    std::cout << "===================" << std::endl;
    std::cout << "===================" << std::endl;
    std::cout << "Getting OMP parameters for Z = " << Z << ", A = " << A_dd << ", energy = " << energy << " MeV"
              << std::endl;
    std::cout << "===================" << std::endl;
    std::cout << "===================" << std::endl;

    // Print OMP parameters for different OMPs
    PhysOMP::Haixia haixia(Z, A_dd, energy * 2);
    std::cout << " ------------------------------------------------" << std::endl;
    std::cout << "Haixia OMP parameters:" << std::endl;
    std::cout << " ------------------------------------------------" << std::endl;
    haixia.GetFrescoOutput(potIdx);
    PhysOMP::Daehnick daehnick(Z, A_dd, energy * 2);
    std::cout << " ------------------------------------------------" << std::endl;
    std::cout << "Daehnick OMP parameters:" << std::endl;
    std::cout << " ------------------------------------------------" << std::endl;
    daehnick.GetFrescoOutput(potIdx);
    PhysOMP::DA1p da1p(Z, A_dd, energy * 2);
    std::cout << " ------------------------------------------------" << std::endl;
    std::cout << "DA1p OMP parameters:" << std::endl;
    std::cout << " ------------------------------------------------" << std::endl;
    da1p.GetFrescoOutput(potIdx);
    PhysOMP::KoningDelaroche koningDelaroche(Z, A_dp, equivalentEnergy);
    std::cout << " ------------------------------------------------" << std::endl;
    std::cout << "Koning-Delaroche OMP parameters:" << std::endl;
    std::cout << " ------------------------------------------------" << std::endl;
    koningDelaroche.GetFrescoOutput(potIdx);
    PhysOMP::CH89 ch89(Z, A_dp, equivalentEnergy);
    std::cout << " ------------------------------------------------" << std::endl;
    std::cout << "CH89 OMP parameters:" << std::endl;
    std::cout << " ------------------------------------------------" << std::endl;
    ch89.GetFrescoOutput(potIdx);
    PhysOMP::Watson watson(Z, A_dp, equivalentEnergy);
    std::cout << " ------------------------------------------------" << std::endl;
    std::cout << "Watson OMP parameters:" << std::endl;
    std::cout << " ------------------------------------------------" << std::endl;
    watson.GetFrescoOutput(potIdx);
}