#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActSilData.h"
#include "ActSilSpecs.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"

#include <fstream>

void DebugF2Multiplicity()
{
    // Get the runs from the data manager reading the config file
    ActRoot::DataManager dataManager {};
    dataManager.ReadDataFile("../../configs/data.conf");

    // For L1 enough with 4 runs (64, 67). For lat sils put at least 10
    // dataManager.SetRuns(60, 70);

    // Get df for the runs
    auto chain {dataManager.GetChain(ActRoot::ModeType::EReadSilMod)};
    auto chain1 {dataManager.GetChain(ActRoot::ModeType::EMerge)};
    chain->AddFriend(chain1.get(), ("MergerData"));
    auto df {ROOT::RDataFrame(*chain)};

    // GATCONFs: 1 -> l0, 2 -> r0, 4 -> f0, 8 -> L1,

    // Get silspecs
    ActPhysics::SilSpecs silspecs {};
    silspecs.ReadFile("../../configs/silspecs.conf");

    auto dfFilterGATCONF = df.Filter(
        [](ActRoot::ModularData& mod)
        {
            if(mod.Get("GATCONF") == 1 || mod.Get("GATCONF") == 2)
                return true;
            return false;
        },
        {"ModularData"});

    // Remember that silData multiplicity is not the same as the mergerData multiplicity

    // Definethe multiplicity of F2
    auto dfMultiplicityF2 = dfFilterGATCONF.Define(
        "F2Multiplicity",
        [&silspecs](ActRoot::SilData& s)
        {
            for(auto layer : s.GetLayers())
            {
                if(layer == "f2")
                {
                    int hits = 0;
                    for(int i = 0; i < s.fSiN[layer].size(); i++)
                    {
                        if(s.fSiE[layer][i] > silspecs.GetLayer("f2").ApplyThreshold(s.fSiN[layer][i], 1.0))
                            hits++;
                    }
                    return hits;
                }
            }
            return -1111;
        },
        {"SilData"});

    auto dfEnergyF2mult = dfMultiplicityF2
                              .Define("EnergyF2mult1",
                                      [](ActRoot::SilData& s, int f2mult)
                                      {
                                          if(f2mult == 1)
                                          {
                                              for(auto layer : s.GetLayers())
                                              {
                                                  if(layer == "f2")
                                                  {
                                                      return s.fSiE[layer][0];
                                                  }
                                              }
                                          }
                                          return -1111.f;
                                      },
                                      {"SilData", "F2Multiplicity"})
                              .Define("EnergyF2_0_MultGreater1",
                                      [](ActRoot::SilData& s, int f2mult)
                                      {
                                          if(f2mult > 1)
                                          {
                                              for(auto layer : s.GetLayers())
                                              {
                                                  if(layer == "f2")
                                                  {
                                                      if(s.fSiE[layer].size() > 0)
                                                          return s.fSiE[layer][0];
                                                      else
                                                          return -1111.f;
                                                  }
                                              }
                                          }
                                          return -1111.f;
                                      },
                                      {"SilData", "F2Multiplicity"})
                              .Define("EnergyF2_1_MultGreater1",
                                      [](ActRoot::SilData& s, int f2mult)
                                      {
                                          if(f2mult > 1)
                                          {
                                              for(auto layer : s.GetLayers())
                                              {
                                                  if(layer == "f2")
                                                  {
                                                      if(s.fSiE[layer].size() > 1)
                                                          return s.fSiE[layer][1];
                                                      else
                                                          return -1111.f;
                                                  }
                                              }
                                          }
                                          return -1111.f;
                                      },
                                      {"SilData", "F2Multiplicity"})
                              .Define("EnergyF2_2_MultGreater1",
                                      [](ActRoot::SilData& s, int f2mult)
                                      {
                                          if(f2mult > 1)
                                          {
                                              for(auto layer : s.GetLayers())
                                              {
                                                  if(layer == "f2")
                                                  {
                                                      if(s.fSiE[layer].size() > 2)
                                                          return s.fSiE[layer][2];
                                                      else
                                                          return -1111.f;
                                                  }
                                              }
                                          }
                                          return -1111.f;
                                      },
                                      {"SilData", "F2Multiplicity"})
                              .Define("EnergyF2_3_MultGreater1",
                                      [](ActRoot::SilData& s, int f2mult)
                                      {
                                          if(f2mult > 1)
                                          {
                                              for(auto layer : s.GetLayers())
                                              {
                                                  if(layer == "f2")
                                                  {
                                                      if(s.fSiE[layer].size() > 3)
                                                          return s.fSiE[layer][3];
                                                      else
                                                          return -1111.f;
                                                  }
                                              }
                                          }
                                          return -1111.f;
                                      },
                                      {"SilData", "F2Multiplicity"});

    // Get events with F2 multiplicity = 1
    auto dfFilterF2Multiplicity1 = dfEnergyF2mult.Filter([](int f2mult) { return f2mult == 1; }, {"F2Multiplicity"});
    auto dfFilterF2MultiplicityGreater1 =
        dfEnergyF2mult.Filter([](int f2mult) { return f2mult > 1; }, {"F2Multiplicity"});

    // Plot the multiplicity of F2
    auto hMultiplicityF2 = dfEnergyF2mult.Histo1D(
        {"hMultiplicityF2", "Multiplicity of F2;Multiplicity;Counts", 5, -0.5, 4.5}, "F2Multiplicity");

    // Plot the multiplicity of F2 for events with F2 multiplicity > 1
    auto hMultiplicityF2Greater1 = dfFilterF2MultiplicityGreater1.Histo1D(
        {"hMultiplicityF2Greater1", "Multiplicity of F2 for events with F2 multiplicity > 1;Multiplicity;Counts", 5,
         -0.5, 4.5},
        "F2Multiplicity");

    // Plot energy for events with F2 multiplicity = 1
    auto hEnergyF2Multiplicity1 = dfFilterF2Multiplicity1.Histo1D(
        {"hEnergyF2Multiplicity1", "Energy of F2 for events with F2 multiplicity = 1;Energy [MeV];Counts", 100, 0, 20},
        "EnergyF2mult1");

    // Get energies for each pad for events with F2 multiplicity > 1
    auto hEnergyF2_0_MultGreater1 = dfFilterF2MultiplicityGreater1.Histo1D(
        {"hEnergyF2_0_MultGreater1", "Energy of F2 pad 0 for events with F2 multiplicity > 1;Energy [MeV];Counts", 100,
         0, 10},
        "EnergyF2_0_MultGreater1");
    auto hEnergyF2_1_MultGreater1 = dfFilterF2MultiplicityGreater1.Histo1D(
        {"hEnergyF2_1_MultGreater1", "Energy of F2 pad 1 for events with F2 multiplicity > 1;Energy [MeV];Counts", 100,
         0, 10},
        "EnergyF2_1_MultGreater1");
    auto hEnergyF2_2_MultGreater1 = dfFilterF2MultiplicityGreater1.Histo1D(
        {"hEnergyF2_2_MultGreater1", "Energy of F2 pad 2 for events with F2 multiplicity > 1;Energy [MeV];Counts", 100,
         0, 10},
        "EnergyF2_2_MultGreater1");
    auto hEnergyF2_3_MultGreater1 = dfFilterF2MultiplicityGreater1.Histo1D(
        {"hEnergyF2_3_MultGreater1", "Energy of F2 pad 3 for events with F2 multiplicity > 1;Energy [MeV];Counts", 100,
         0, 10},
        "EnergyF2_3_MultGreater1");

    auto* cMultiplicityF2 = new TCanvas("cMultiplicityF2", "Multiplicity of F2", 800, 600);
    cMultiplicityF2->Divide(2, 1);
    cMultiplicityF2->cd(1);
    hMultiplicityF2->DrawClone();
    cMultiplicityF2->cd(2);
    hMultiplicityF2Greater1->DrawClone();

    auto* cEnergyF2Multiplicity1 =
        new TCanvas("cEnergyF2Multiplicity1", "Energy of F2 for events with F2 multiplicity = 1", 800, 600);
    hEnergyF2Multiplicity1->DrawClone();

    auto* cEnergyF2_MultGreater1 =
        new TCanvas("cEnergyF2_0_MultGreater1", "Energy of F2 all pads for events with F2 multiplicity > 1", 800, 600);
    cEnergyF2_MultGreater1->Divide(2, 2);
    cEnergyF2_MultGreater1->cd(1);
    hEnergyF2_0_MultGreater1->DrawClone();
    cEnergyF2_MultGreater1->cd(2);
    hEnergyF2_1_MultGreater1->DrawClone();
    cEnergyF2_MultGreater1->cd(3);
    hEnergyF2_2_MultGreater1->DrawClone();
    cEnergyF2_MultGreater1->cd(4);
    hEnergyF2_3_MultGreater1->DrawClone();
}