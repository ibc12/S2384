#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include "./GetSigma.cxx"
#include "./do_simuL1.cxx"
#include "./plotter.cxx"

void runnerL1(TString what = "simu", bool inspect = false)
{
    // Neutron and Proton phase space
    int neutronPS {0}; // number of neutrons in final state
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
            // Exs = {0, 0.981};
            // Exs = {4.100, 5.400, 6.100, 6.530, 7.100};
            // Exs = {2.255};
            // Exs = {0.};
            // Exs = {0, 0.981, 2.255, 3.210, 4.100, 5.400, 6.100, 6.530, 7.100, 8.000, 9.000, 10.000};
            Exs = {0.000, 2.700, 4.600, 8.100, 8.600};
    }

    else if(neutronPS == 0 && protonPS == 0 && target == "2H" && light == "2H") // Elastic and Inelastic scattering
    {
        if(beam == "11Li")
            // Exs = {0, 1.266, 2.474};
            Exs = {0};
        else if(beam == "7Li")
            // Exs = {0, 0.477};
            Exs = {0}; // Only g.s.
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

    // Output directory and file naming (single source of truth)
    const TString outDir {TString::Format("./Outputs/%s", beam.c_str())};
    // suffix = "" -> final file; suffix = "_part<i>" -> partial file of job i
    auto outName {[&](double ex, const TString& suffix = "")
                  {
                      return TString::Format("%s/%s_%s_TRIUMF_Eex_%.3f_nPS_%d_pPS_%d_L1%s.root", outDir.Data(),
                                             target.c_str(), light.c_str(), ex, neutronPS, protonPS, suffix.Data());
                  }};

    // Run simu or plot
    ROOT::EnableThreadSafety();
    if(what.Contains("simu"))
    {
        gSystem->mkdir(outDir.Data(), true); // recursive, does nothing if it already exists

        if(inspect)
        {
            do_simuL1(beam, target, light, heavy, neutronPS, protonPS, Tbeam, Exs.front(), inspect,
                      outName(Exs.front()).Data());
        }
        else
        {
            // Number of parts in which each Ex is split (granularity of the work)
            const int nParts {isPS ? 8 : 2};
            // Maximum number of simultaneous ROOT processes (global, across all Ex)
            const int maxJobs {12};

            // Task queue: one task per (Ex, part)
            struct Task
            {
                size_t iEx;
                int part;
            };
            std::vector<Task> tasks;
            tasks.reserve(Exs.size() * nParts);
            for(size_t iEx = 0; iEx < Exs.size(); iEx++)
                for(int i = 1; i <= nParts; i++)
                    tasks.push_back({iEx, i});

            std::vector<int> remaining(Exs.size(), nParts); // pending parts per Ex
            std::mutex mtx;                                 // protects remaining and std::cout
            std::atomic<size_t> next {0};                   // index of the next task to take

            auto poolWorker {
                [&]()
                {
                    while(true)
                    {
                        const size_t t {next++};
                        if(t >= tasks.size())
                            break;
                        const size_t iEx {tasks[t].iEx};
                        const double ex {Exs[iEx]};
                        const TString partFile {outName(ex, TString::Format("_part%d", tasks[t].part))};
                        const TString cmd {
                            TString::Format("root -l -b -q -x "
                                            "'do_simuL1.cxx(\"%s\",\"%s\",\"%s\",\"%s\",%d,%d,%f,%f,%d,\"%s\")'",
                                            beam.c_str(), target.c_str(), light.c_str(), heavy.c_str(), neutronPS,
                                            protonPS, Tbeam, ex, inspect, partFile.Data())};
                        std::system(cmd.Data());

                        // Did I finish the last part of this Ex?
                        bool last {false};
                        {
                            std::lock_guard<std::mutex> lock {mtx};
                            last = (--remaining[iEx] == 0);
                        }
                        if(!last)
                            continue;

                        // Merge and clean this Ex without waiting for the others
                        const TString haddout {outName(ex)};
                        TString haddlist {};
                        for(int i = 1; i <= nParts; i++)
                            haddlist += outName(ex, TString::Format("_part%d", i)) + " ";
                        {
                            std::lock_guard<std::mutex> lock {mtx};
                            std::cout << "Output file: " << haddout << '\n';
                            std::cout << "Input files: " << haddlist << '\n';
                        }
                        const int rc {
                            std::system(TString::Format("hadd -f %s %s", haddout.Data(), haddlist.Data()).Data())};
                        if(rc == 0)
                            std::system(TString::Format("rm -f %s", haddlist.Data()).Data());
                        else
                        {
                            std::lock_guard<std::mutex> lock {mtx};
                            std::cerr << "hadd failed for Ex = " << ex << ". Partial files kept for inspection\n";
                        }
                    }
                }};

            const int nThreads {std::min<int>(maxJobs, static_cast<int>(tasks.size()))};
            std::vector<std::thread> pool;
            pool.reserve(nThreads);
            for(int i = 0; i < nThreads; i++)
                pool.emplace_back(poolWorker);
            for(auto& th : pool)
                th.join();
        }
    }
    else if(what.Contains("plot"))
    {
        Plotter(Exs, beam, target, light, Tbeam, neutronPS, protonPS);
    }
    else if(what.Contains("sigma"))
    {
        GetSigma(beam, light, target, Exs, neutronPS, protonPS, true);
    }
    else
    {
        std::cout << "That method was not implemented yet" << '\n';
    }
}