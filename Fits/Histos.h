#ifndef Histos_h
#define Histos_h
#include "ROOT/RDF/HistoModels.hxx"

#include "TString.h"
namespace S2384Fit
{
// (d,d) settings
const ROOT::RDF::TH1DModel Exdd {
    "hEx", TString::Format("^{11}Li(d,d);E_{x} [MeV];Counts / %.0f keV", (25. - (-5.)) / 200 * 1000), 200, -5, 25};
// (d,t) settings
const ROOT::RDF::TH1DModel Exdt {
    "hEx", TString::Format("^{11}Li(d,t);E_{x} [MeV];Counts / %.0f keV", (25. - (-5.)) / 100 * 1000), 100, -5, 25};
} // namespace S2384Fit
#endif
