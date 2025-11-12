#include "TString.h"

#include <string>
#include <unordered_map>
#include <vector>

#include "./do_simu.cxx"

void runner(TString what = "simu", bool inspect = true)
{
    // Beam energy
    double Tbeam {7 * 7.558}; // MeV
    // Neutron and Proton phase space
    int neutronPS {0}; // number of neutrons in final state
    int protonPS {0};  // number of protons in final state
    // Particles
    std::string beam {"7Li"};
    std::string target {"2H"};
    std::string light {"1H"};
    std::string heavy {"8Li"};
    // Vector with Exs
    std::vector<double> Exs;
    if(neutronPS == 0 && protonPS == 0 && target == "2H" && light == "1H") // Transfer dp
    {
        if(beam == "11Li")
            Exs = {0, 0.130, 0.435};
        else if(beam == "7Li")
            Exs = {0, 0.981, 2.255};
    }
        
    else if(neutronPS == 0 && protonPS == 0 && target == "2H" && light == "2H") // Elastic and Inelastic scattering
    {
        if(beam == "11Li")
            Exs = {0, 1.266, 2.474};
        else if(beam == "7Li")
            Exs = {0, 0.477};
    }
    else if(target == "2H" && light == "3H") // dt (only g.s)
        Exs = {0};
    else if(heavy == "9He" || heavy == "10He") // d3He,alpha (only g.s)
        Exs = {0};
    else if(neutronPS == 2 && protonPS == 0 && target == "2H" && light == "2H")
        Exs = {(1.26642 + 0.36928) / 2}; // half value between first excited state and the S_2n
    else if(neutronPS > 0 && protonPS == 0)
        Exs = {0}; // only gs for n phase space
    else if(neutronPS == 0 && protonPS > 0)
        Exs = {0};
    else
        throw std::runtime_error("No confs with neutronPS and protonPS enabled at the same time");

    // Run simu or plot
    if(what.Contains("simu"))
    {
        for(const auto& ex : Exs)
        {
            do_all_simus(beam, target, light, heavy, neutronPS, protonPS ,Tbeam, ex, inspect);
            // auto str {TString::Format("root -l -b -x -q \'triumf.cxx(\"%s\",\"%s\",\"%s\",%f,%f,%d)\'", beam.c_str(),
            //                           target.c_str(), light.c_str(), Tbeam, ex, inspect)};
            if(inspect)
                break; // inspect: to debug simulation
        }
    }
    else
    {
        std::cout << "That method was not implemented yet" << '\n';
    }
}
