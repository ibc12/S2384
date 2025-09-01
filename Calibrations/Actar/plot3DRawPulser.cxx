#include "TChain.h"
#include "TString.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TColor.h"

#include "ActCalibrationManager.h"
#include "ActTPCDetector.h"

void FillHistogram(ActRoot::CalibrationManager *calman, ActRoot::TPCParameters *pars, MEventReduced *evt, TH2D *h)
{

    // iterate over hits
    for (int it = 0, size = evt->CoboAsad.size(); it < size; it++)
    {
        int co = evt->CoboAsad[it].globalchannelid >> 11;
        int as = (evt->CoboAsad[it].globalchannelid - (co << 11)) >> 9;
        int ag = (evt->CoboAsad[it].globalchannelid - (co << 11) - (as << 9)) >> 7;
        int ch = evt->CoboAsad[it].globalchannelid - (co << 11) - (as << 9) - (ag << 7);
        int where = co * pars->GetNBASAD() * pars->GetNBAGET() * pars->GetNBCHANNEL() +
                    as * pars->GetNBAGET() * pars->GetNBCHANNEL() + ag * pars->GetNBCHANNEL() + ch;

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
                        // Fill histogram

                        if (Qiaux >= 1500) // in this experiment all runs have a baseline that we eliminate
                        {
                            h->Fill(where, Qiaux);
                        }
                    }
                }
            }
        }
    }
}

void plot3DRawPulser()
{

    // Get the data into TChain
    auto chain{new TChain("ACTAR_TTree")};
    std::vector<int> runs{85}; // Put whichever runs needed (this was for e796)
    for (const auto &run : runs)
    {
        chain->Add(TString::Format("../../RootFiles/Raw/Pulser_ok/Tree_Run_%04d_Merged.root", run));
    }

    // Set the parameters
    ActRoot::TPCParameters tpc{"Actar"};
    // Calibration manager
    ActRoot::CalibrationManager calman{};
    calman.ReadLookUpTable("../Actar/LT.txt"); // convert LT file in a string line

    // Create histogram
    auto *h{new TH2D{"h", "pads;Channel number;Calibrated charge [u.a.]", 17408, 0, 17408, 800, 0, 5000}};

    // Set MEventReduced
    MEventReduced *evt{new MEventReduced};
    chain->SetBranchAddress("data", &evt); // get the column from the chain that we gonna get data from

    for (long int entry = 0, maxEntry = chain->GetEntries(); entry < maxEntry; entry++)
    {
        std::cout << "\r"
                  << "At entry : " << entry << std::flush;
        chain->GetEntry(entry); // get  the data from the chain and write it to evt variable
        FillHistogram(&calman, &tpc, evt, h);
    }

    // Plot
    gStyle->SetOptStat(0);
    auto *c0{new TCanvas{"c0", "Gain matching canvas"}};
    h->GetXaxis()->SetLabelSize(0.03);   // tamaño de los números en el eje X
    h->GetYaxis()->SetLabelSize(0.03);   // tamaño de los números en el eje Y

    h->GetXaxis()->SetTitleSize(0.04);   // tamaño del título del eje X
    h->GetYaxis()->SetTitleSize(0.04);
    h->GetYaxis()->SetTitleOffset(1);  // separa el título del eje Y
    h->GetXaxis()->SetTitleOffset(0.7);
    h->SetFillColor(kAzure+1);  // pick any ROOT color
    h->Draw("BOX");
    auto *c1{new TCanvas{"c1", "Gain matching canvas 2"}};
    h->DrawClone("cloz");

    // Lets get the mean values of y for each channel
    std::vector<std::pair<int,double>> channelMeans; 

    int nChannels = h->GetNbinsX();
    std::cout << "nX: " << nChannels << std::endl;
    int nY = h->GetNbinsY();
    std::cout << "nY: " << nY << std::endl;

    for (int ix = 1; ix <= nChannels; ix++)  // bins empiezan en 1
    {
        double sumW = 0;   // suma de contenidos
        double sumYW = 0;  // suma de y * contenido

        for (int iy = 1; iy <= nY; iy++)
        {
            double content = h->GetBinContent(ix, iy);
            if (content > 0)
            {
                double yval = h->GetYaxis()->GetBinCenter(iy);
                sumYW += yval * content;
                sumW  += content;
            }
        }

        if (sumW > 0)
        {
            double meanY = sumYW / sumW;
            channelMeans.emplace_back(ix, meanY);
        }
        else{
            channelMeans.emplace_back(ix, 0.0); // No data for this channel
        }
    }

    int nx = 128, ny = 128;
    TH2F *h2 = new TH2F("h2", "Mapa XY con amplitud en Z", nx, 0, nx, ny, 0, ny);
    h2->GetXaxis()->SetTitle("X position [pads]");
    h2->GetYaxis()->SetTitle("Y position [pads]");
    h2->GetZaxis()->SetTitle("Mean charge [u.a.]");

    // Opcional: ajustar tamaños y offsets
    h2->GetXaxis()->SetTitleSize(0.04);
    h2->GetYaxis()->SetTitleSize(0.04);
    h2->GetZaxis()->SetTitleSize(0.04);

    h2->GetXaxis()->SetTitleOffset(1.4);
    h2->GetYaxis()->SetTitleOffset(1.4);
    h2->GetZaxis()->SetTitleOffset(1.4);

    // Do lego plot for the run of max charge deposition
    for(int i = 0; i < nChannels; i++)
    {
        auto x = calman.ApplyLookUp(i, 4);
        if(x == -1) 
            continue; // FPN channels
        auto y = calman.ApplyLookUp(i, 5);
        int binx = h2->GetXaxis()->FindBin(x);
        int biny = h2->GetYaxis()->FindBin(y);
        h2->SetBinContent(binx, biny, channelMeans[i].second);
        //h2->SetBinContent(x, y, channelMeans[i].second);
        // if(i % 100 == 0)
        //     std::cout << "Channel: " << i << " X: " << x << " Y: " << y << " Mean charge: " << channelMeans[i].second << std::endl;
    }
    // Dibujar con barras tipo lego 
    h2->SetMinimum(1800);
    h2->SetMaximum(2800);
    gStyle->SetPalette(kGreyScale);
    h2->SetLineColor(kBlack);
    h2->SetLineWidth(1);
    std::cout<<"n bins x"<<h2->GetNbinsX()<<std::endl;
    h2->DrawClone("lego2z");
}