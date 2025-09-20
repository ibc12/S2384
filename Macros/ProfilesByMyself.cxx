#include "ActVoxel.h"
#include "ActCluster.h"

#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"
#include "TMath.h"
#include "TF1.h"

void ScalePoint(ROOT::Math::XYZPointF &point, float xy, float z, bool addOffset = false) 	
{
    if(addOffset) // when converting a bin point to physical units which wasnt already corrected
        point += ROOT::Math::XYZVector {0.5, 0.5, 0.5};
    point.SetX(point.X() * xy);
    point.SetY(point.Y() * xy);
    point.SetZ(point.Z() * z);
}

TH1D* ProfilesByActRoot(const ActRoot::TPCData &d, const ActRoot::MergerData &m)
{
    auto h = new TH1D("h", "Perfil de carga", 180, -5, 265);

        auto cluster = d.fClusters[m.fLightIdx];
        auto voxels = cluster.GetRefToVoxels();

        // Punto de proyección del voxel más frontal como ActRoot
        auto front {voxels.front().GetPosition()};
        auto r0 {cluster.GetLine().ProjectionPointOnLine(front)};
        ScalePoint(r0, 2., 2.84032, true);
        
        auto line = cluster.GetLine();
        line.Scale(2., 2.84032); // Escalar a mm
        line.AlignUsingPoint(r0, true);
        auto dir = line.GetDirection().Unit();
        // Divisiones por voxel
        float div = 1.f/3;
        for (auto &v : voxels)
        {
            auto pos = v.GetPosition();
            double q = v.GetCharge();

            for (int ix=-1; ix<=1; ix++)
                for (int iy=-1; iy<=1; iy++)
                    for (int iz=-1; iz<=1; iz++)
                    {
                        ROOT::Math::XYZPointF bin {(pos.X() + 0.5f) + ix * div, (pos.Y() + 0.5f) + iy * div,
                                    (pos.Z() + 0.5f) + iz * div};
                        ScalePoint(bin, 2., 2.84032);
                        // Convert to physical units
                        // Project it on line
                        auto proj {line.ProjectionPointOnLine(bin)};
                        auto position {(proj - r0).R()};
                        h->Fill(position, q / 27.0); // 27 mini-voxels
                    }
        }
    return h;
}

TH1D* SmoothHistogramPreserveIntegral(TH1D* h, int nBinsKernel = 5, double sigma = 1.0)
{
    if (!h) return nullptr;

    // Store original integral
    double integral = h->Integral();

    // Create Gaussian kernel
    std::vector<double> kernel(nBinsKernel);
    int half = nBinsKernel / 2;
    double sumKernel = 0;
    for (int i = 0; i < nBinsKernel; ++i) {
        double x = i - half;
        kernel[i] = TMath::Exp(-0.5 * x * x / (sigma * sigma));
        sumKernel += kernel[i];
    }
    for (auto &k : kernel) k /= sumKernel; // normalize kernel

    // Convolve histogram
    TH1D* hSmooth = (TH1D*)h->Clone(Form("%s_smooth", h->GetName()));
    hSmooth->Reset();
    int nBins = h->GetNbinsX();
    for (int i = 1; i <= nBins; ++i) {
        double sum = 0;
        for (int j = 0; j < nBinsKernel; ++j) {
            int bin = i + j - half;
            if (bin >= 1 && bin <= nBins)
                sum += h->GetBinContent(bin) * kernel[j];
        }
        hSmooth->SetBinContent(i, sum);
    }

    // Preserve integral
    hSmooth->Scale(integral / hSmooth->Integral());

    return hSmooth;
}


void ProfilesByMyself()
{
    std::string dataconf{"./../configs/data.conf"};
    ActRoot::DataManager dataman{dataconf, ActRoot::ModeType::EFilter};
    auto chain{dataman.GetChain()};
    auto chain2{dataman.GetChain(ActRoot::ModeType::EMerge)};
    chain->AddFriend(chain2.get());

    // RDataFrame
    //ROOT::EnableImplicitMT();
    ROOT::RDataFrame dforigin{*chain};

    auto df = dforigin.Filter(
        [](ActRoot::MergerData &m)
        {
            // Filter condition here
            if (m.fRun == 19 && m.fEntry == 519)
                return true;
            else
                return false;
        },
        {"MergerData"});

    auto df2 = df.Define("ChargeProfile", [](const ActRoot::TPCData &d, const ActRoot::MergerData &m)
    {
        auto h = new TH1D("h", "Perfil de carga", 180*2, -265, 265);

        auto cluster = d.fClusters[m.fLightIdx];
        auto line = cluster.GetLine();
        line.Scale(2., 2.84032); // Escalar a mm
        auto dir = line.GetDirection().Unit();
        cluster.SortAlongDir(dir);
        auto voxels = cluster.GetRefToVoxels();
        

        // Centro de gravedad cargado
        // double totalCharge = 0;
        // double x=0, y=0, z=0;
        // for (auto &v : voxels) {
        //     auto p = v.GetPosition();
        //     //ScalePoint(p, 2., 2.84032);
        //     double q = v.GetCharge();
        //     x += p.X() * q;
        //     y += p.Y() * q;
        //     z += p.Z() * q;
        //     totalCharge += q;
        // }
        // ROOT::Math::XYZPointF r0(x/totalCharge, y/totalCharge, z/totalCharge);
        

        // Punto de proyección del voxel más frontal como ActRoot
        auto front {voxels.front().GetPosition()};
        auto r0 {cluster.GetLine().ProjectionPointOnLine(front)};
        std::cout << "Front at: (" << front.X() << ", " << front.Y() << ", " << front.Z() << ")\n";
        ScalePoint(r0, 2., 2.84032, true);
        std::cout << "Centroid at: (" << r0.X() << ", " << r0.Y() << ", " << r0.Z() << ")\n";

        line.AlignUsingPoint(r0, true);
        auto r0OnLine {line.MoveToX(r0.X())};
        std::cout<< "r0 on line: (" << r0OnLine.X() << ", " << r0OnLine.Y() << ", " << r0OnLine.Z() << ")\n";
        // Divisiones por voxel
        float div = 1.f/3;
        for (auto &v : voxels)
        {
            auto pos = v.GetPosition();
            double q = v.GetCharge();

            for (int ix=-1; ix<=1; ix++)
                for (int iy=-1; iy<=1; iy++)
                    for (int iz=-1; iz<=1; iz++)
                    {
                        // Miña maneira
                        //ROOT::Math::XYZPointF miniPoint(pos.X() + ix*div,
                        //                                pos.Y() + iy*div,
                        //                                pos.Z() + iz*div);
                        //// proyectar sobre el cluster
                        //auto r = miniPoint - r0;
                        // double position = r.Dot(dir);
                        ROOT::Math::XYZPointF bin {(pos.X() + 0.5f) + ix * div, (pos.Y() + 0.5f) + iy * div,
                                    (pos.Z() + 0.5f) + iz * div};
                        ScalePoint(bin, 2., 2.84032);
                        // Convert to physical units
                        // Project it on line
                        auto proj {line.ProjectionPointOnLine(bin)};
                        auto position {(proj - r0).Dot(dir)}; // can bealso done with .R() but if voxel not first in direction "dir" it gives accumulation for low distances
                        h->Fill(position, q / 27.0); // 27 mini-voxels
                    }
        }
        // h->Rebin(2);
        TH1D* hSmooth = SmoothHistogramPreserveIntegral(h, 5, 1.0);
        //delete h; // optional: free original histogram if not needed
        return h;
    }, {"TPCData", "MergerData"});

    auto profiles = df2.Take<TH1D*>("ChargeProfile");
    profiles->at(0)->DrawClone("hist");
}