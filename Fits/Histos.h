#ifndef Histos_h
#define Histos_h
#include "ROOT/RDF/HistoModels.hxx"

#include "TString.h"
namespace S2384Fit
{
// (d,p) settings 7Li
const ROOT::RDF::TH1DModel Exdp_7Li {
    "hEx", TString::Format("^{7}Li(d,p);E_{x} [MeV];Counts / %.0f keV", (25. - (-5.)) / 300 * 1000), 300, -5, 25};
// (d,d) settings 7Li
const ROOT::RDF::TH1DModel Exdd_7Li {
    "hEx", TString::Format("^{7}Li(d,d);E_{x} [MeV];Counts / %.0f keV", (10. - (-5.)) / 200 * 1000), 200, -5, 10};
// (d,t) settings 7Li
const ROOT::RDF::TH1DModel Exdt_7Li {
    "hEx", TString::Format("^{7}Li(d,t);E_{x} [MeV];Counts / %.0f keV", (25. - (-5.)) / 100 * 1000), 100, -5, 25};
// (d,p) settings 11Li
const ROOT::RDF::TH1DModel Exdp {
    "hEx", TString::Format("^{11}Li(d,p);E_{x} [MeV];Counts / %.0f keV", (25. - (-5.)) / 200 * 1000), 200, -5, 25};
// (d,d) settings 11Li
const ROOT::RDF::TH1DModel Exdd {
    "hEx", TString::Format("^{11}Li(d,d);E_{x} [MeV];Counts / %.0f keV", (25. - (-5.)) / 200 * 1000), 200, -5, 25};
// (d,t) settings 11Li
const ROOT::RDF::TH1DModel Exdt {
    "hEx", TString::Format("^{11}Li(d,t);E_{x} [MeV];Counts / %.0f keV", (25. - (-5.)) / 100 * 1000), 100, -5, 25};
} // namespace S2384Fit
#endif
