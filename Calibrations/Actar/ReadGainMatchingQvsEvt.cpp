#include "TChain.h"
#include "TString.h"
#include "TH2D.h"
#include "TCanvas.h"

#include "ActCalibrationManager.h"
#include "ActTPCDetector.h"

void FillHistogram(ActRoot::CalibrationManager *calman, ActRoot::TPCParameters *pars,
                   MEventReduced *evt, TH2D *h, bool isMatched,
                   TH2D *hQvsEntryNumber = nullptr, long int entry = -1)
{
    for (int it = 0, size = evt->CoboAsad.size(); it < size; it++)
    {
        int co = evt->CoboAsad[it].globalchannelid >> 11;
        int as = (evt->CoboAsad[it].globalchannelid - (co << 11)) >> 9;
        int ag = (evt->CoboAsad[it].globalchannelid - (co << 11) - (as << 9)) >> 7;
        int ch = evt->CoboAsad[it].globalchannelid - (co << 11) - (as << 9) - (ag << 7);
        int where = co * pars->GetNBASAD() * pars->GetNBAGET() * pars->GetNBCHANNEL() +
                    as * pars->GetNBAGET() * pars->GetNBCHANNEL() +
                    ag * pars->GetNBCHANNEL() + ch;

        if ((co != 31) && (co != 16))
        {
            auto xval{calman->ApplyLookUp(where, 4)};
            auto yval{calman->ApplyLookUp(where, 5)};
            for (int hit = 0, otherSize = evt->CoboAsad[it].peakheight.size(); hit < otherSize; hit++)
            {
                if ((yval != -1) && (xval != -1))
                {
                    double z_position{evt->CoboAsad[it].peaktime[hit]};
                    if (z_position > 0.)
                    {
                        auto Qiaux{evt->CoboAsad[it].peakheight[hit]};
                        if (isMatched)
                            Qiaux = calman->ApplyPadAlignment(where, Qiaux);

                        if (Qiaux >= 0)
                        {
                            h->Fill(where, Qiaux);
                            if (hQvsEntryNumber && entry >= 0)
                                hQvsEntryNumber->Fill(entry, Qiaux);
                        }
                    }
                }
            }
        }
    }
}

void ReadGainMatchingQvsEvt(bool isMatched = true)
{
    // Load TChain
    auto chain{new TChain("ACTAR_TTree")};
    std::vector<int> runs{85}; // Update with your run list
    for (const auto &run : runs)
    {
        chain->Add(TString::Format("../../RootFiles/Raw/Tree_Run_%04d_Merged.root", run));
    }

    // Set TPC parameters and calibration manager
    ActRoot::TPCParameters tpc{"Actar"};
    ActRoot::CalibrationManager calman{};
    calman.ReadLookUpTable("../Actar/LT.txt");
    if (isMatched)
        calman.ReadPadAlign("./Outputs/gain_matching_s2384_v0.dat");

    // Histograms
    auto *h{new TH2D{"h", "pads;Channel;Q", 17408, 0, 17408, 800, 0, 5000}};
    auto *hQvsEntryNumber{new TH2D{"hQvsEntryNumber", "Q vs Entry Number;Entry Number;Q", 1000, -2, 1000, 2500, -2, 2500}};

    // Set event pointer
    MEventReduced *evt{new MEventReduced};
    chain->SetBranchAddress("data", &evt);

    // Event loop
    long int maxEntry = chain->GetEntries();
    for (long int entry = 0; entry < maxEntry; ++entry)
    {
        std::cout << "\rProcessing entry: " << entry << "/" << maxEntry << std::flush;
        chain->GetEntry(entry);
        FillHistogram(&calman, &tpc, evt, h, isMatched, hQvsEntryNumber, entry);
    }
    std::cout << "\nDone.\n";

    // Plot gain matching histogram
    auto *c0{new TCanvas{"c0", "Gain matching canvas", 1200, 800}};
    h->DrawClone("colz");

    // Plot Q vs Entry Number histogram
    auto *c1{new TCanvas{"c1", "Q vs Entry Number", 1200, 800}};
    hQvsEntryNumber->DrawClone("colz");

    // Save histograms if desired
    if (!isMatched)
        h->SaveAs("./Inputs/gain_85.root");

    hQvsEntryNumber->SaveAs("./Outputs/Q_vs_EntryNumber.root");
}
