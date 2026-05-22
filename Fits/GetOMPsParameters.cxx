#include "PhysOMP.h"

#include <iostream>

void GetOMPsParameters()
{
    int Z = 3;               // Lithium
    int A = 8;               // Lithium-11
    double energy {};
    if(Z == 3)
    {
        energy = 7.26557; // 7Li beam time
    }
    else if(Z == 3)
    {
        energy = 7.37188; // 11Li beam time
    }
    std::cout << "===================" << std::endl;
    std::cout << "===================" << std::endl;
    std::cout << "Getting OMP parameters for Z = " << Z << ", A = " << A << ", energy = " << energy << " MeV" << std::endl;
    std::cout << "===================" << std::endl;
    std::cout << "===================" << std::endl;

    // Print OMP parameters for different OMPs
    PhysOMP::Haixia haixia(Z, A, energy);
    std::cout << " ------------------------------------------------" << std::endl;
    std::cout << "Haixia OMP parameters:" << std::endl;
    std::cout << " ------------------------------------------------" << std::endl;
    haixia.Print();
    PhysOMP::Daehnick daehnick(Z, A, energy);
    std::cout << " ------------------------------------------------" << std::endl;
    std::cout << "Daehnick OMP parameters:" << std::endl;
    std::cout << " ------------------------------------------------" << std::endl;
    daehnick.Print();
    PhysOMP::DA1p da1p(Z, A, energy);
    std::cout << " ------------------------------------------------" << std::endl;
    std::cout << "DA1p OMP parameters:" << std::endl;
    std::cout << " ------------------------------------------------" << std::endl;
    da1p.Print();
    PhysOMP::KoningDelaroche koningDelaroche(Z, A, energy);
    std::cout << " ------------------------------------------------" << std::endl;
    std::cout << "Koning-Delaroche OMP parameters:" << std::endl;
    std::cout << " ------------------------------------------------" << std::endl;
    koningDelaroche.Print();
    PhysOMP::CH89 ch89(Z, A, energy);
    std::cout << " ------------------------------------------------" << std::endl;
    std::cout << "CH89 OMP parameters:" << std::endl;
    std::cout << " ------------------------------------------------" << std::endl;
    ch89.Print();
}