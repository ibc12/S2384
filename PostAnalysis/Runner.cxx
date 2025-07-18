#include "ActColors.h"

#include "TROOT.h"
#include "TString.h"

#include <iostream>
#include <string>

void Runner(TString what = "")
{
    std::string beam {"11Li"};
    std::string target {"d"};
    std::string light {"p"};

    std::cout << BOLDGREEN << "···· Runner ····" << '\n';
    std::cout << "-> Beam   : " << beam << '\n';
    std::cout << "-> Target : " << target << '\n';
    std::cout << "-> Light  : " << light << '\n';
    std::cout << "-> What   : " << what << '\n';
    std::cout << "······························" << RESET << '\n';

    auto args {TString::Format("(\"%s\", \"%s\", \"%s\")", beam.c_str(), target.c_str(), light.c_str())};
    TString path {"./Pipes/"};
    TString func {};
    TString ext {".cxx"};

    // PID
    if(what.Contains("1"))
    {
        func = "Pipe1_PID";
        gROOT->LoadMacro(path + func + ext);
        gROOT->ProcessLine(func + args);
    }
}
