#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"

#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include "./do_simu.cxx"
#include "./plotter.cxx"

void runner(TString what = "simu", bool inspect = true)
{
    // Neutron and Proton phase space
    int neutronPS {1}; // number of neutrons in final state
    int protonPS {0};  // number of protons in final state
    bool isPS {neutronPS > 0 || protonPS > 0};
    // Particles
    std::string beam {"7Li"};
    std::string target {"2H"};
    std::string light {"2H"};
    std::string heavy {"7Li"};
    // Beam energy
    double Tbeam {};
    if(beam == "7Li")
        Tbeam = 7 * 7.558; // MeV
    else if(beam == "11Li")
        Tbeam = 11 * 7.558; // MeV
    // Vector with Exs
    std::vector<double> Exs;
    if(neutronPS == 0 && protonPS == 0 && target == "2H" && light == "1H") // Transfer dp
    {
        if(beam == "11Li")
            Exs = {0, 0.130, 0.435};
        else if(beam == "7Li")
            // Exs = {0,  0.981, 2.255, 3.210};
            // Exs = {5.400, 6.100, 6.530, 7.100};
            Exs = {0, 0.981, 2.255, 3.210, 5.400, 6.100, 6.530, 7.100};
    }

    else if(neutronPS == 0 && protonPS == 0 && target == "2H" && light == "2H") // Elastic and Inelastic scattering
    {
        if(beam == "11Li")
            Exs = {0, 1.266, 2.474};
        else if(beam == "7Li")
            Exs = {0, 0.477};
    }
    else if(target == "2H" && light == "3H") // dt (only g.s)
        Exs = {0, 1, 2, 3};
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
    ROOT::EnableThreadSafety();
    std::vector<std::thread> threads;
    auto worker {[](TString str) { return gSystem->Exec(str.Data()); }};
    if(what.Contains("simu"))
    {
        if(inspect)
        {
            do_simu(beam, target, light, heavy, neutronPS, protonPS, Tbeam, Exs.front(), inspect);
        }
        else
        {
            TString haddlist {};
            TString haddout {};
            if(isPS)
            {
                haddout = TString::Format("./Outputs/%s/%s_%s_TRIUMF_Eex_%.3f_nPS_%d_pPS_%d.root", beam.c_str(),
                                          target.c_str(), light.c_str(), Exs.front(), neutronPS, protonPS);
                // List of files generated per thread
                // Number of threads = 6
                int nthreads {6};
                // 1e8 events each
                for(int i = 1; i <= nthreads; i++)
                {
                    auto str {TString::Format(
                        "root -l -b -q -x 'do_simu.cxx(\"%s\",\"%s\",\"%s\",\"%s\",%d,%d,%f,%f,%d,%d)\'",
                        beam.c_str(), target.c_str(), light.c_str(), heavy.c_str(), neutronPS, protonPS, Tbeam,
                        Exs.front(), inspect, i)};
                    haddlist += TString::Format("./Outputs/%s/%s_%s_TRIUMF_Eex_%.3f_nPS_%d_pPS_%d_%s.root",
                                                beam.c_str(), target.c_str(), light.c_str(), Exs.front(), neutronPS,
                                                protonPS, std::to_string(i).c_str()) +
                                " ";
                    threads.emplace_back(worker, str);
                }
            }
            else
            {
                for(const auto& ex : Exs)
                {
                    auto str = TString::Format("root -l -b -q -x 'do_simu.cxx(\"%s\",\"%s\",\"%s\",\"%s\",%d,%d,%f,%f,%d)'",
                                               beam.c_str(), target.c_str(), light.c_str(), heavy.c_str(), neutronPS,
                                               protonPS, Tbeam, ex, inspect);
                    threads.emplace_back(std::thread {worker, str});
                }
            }
            // Join threads for PS simu
            for(auto& th : threads)
                th.join();
            // Now hadd if PS
            if(isPS)
            {
                std::cout << "Output file: " << haddout << '\n';
                std::cout << "Input files: " << haddlist << '\n';
                // Merge
                gSystem->Exec(TString::Format("hadd -f %s %s", haddout.Data(), haddlist.Data()).Data());
                // Remove
                gSystem->Exec(TString::Format("rm %s", haddlist.Data()).Data());
            }
        }
    }
    else if(what.Contains("plot"))
    {
        Plotter(Exs, beam, target, light, Tbeam, neutronPS, protonPS);
    }
    else
    {
        std::cout << "That method was not implemented yet" << '\n';
    }
}
