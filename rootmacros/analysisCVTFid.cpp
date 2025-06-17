#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TCanvas.h>
#include <TFile.h>
#include <TLegend.h>
#include <TMath.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TLatex.h>
#include <iostream>
#include <map>
#include <vector>
#include <tuple>

using namespace ROOT::VecOps;

int thetaRegionIndex(float thetaRad, const std::vector<float> &thetaCuts) {
    float deg = thetaRad * 180.0 / TMath::Pi();
    for (size_t i = 0; i < thetaCuts.size(); ++i) {
        if (deg <= thetaCuts[i]) return i;
    }
    return thetaCuts.size();
}

void DrawCVTChi2ndf_Optimized(const int &selectedPid, const int &selecteddetector,
                               const double ymin, const double ymax, const int nbinsx,
                               const std::vector<int> &layers,
                               const std::vector<float> &xmins,
                               const std::vector<float> &xmaxs,
                               const std::vector<float> &thetaCuts,
                               const std::string &filename, const std::string &treename) {

    ROOT::RDataFrame df(treename, filename);
    gStyle->SetOptStat(0);

    std::vector<std::tuple<int, float, float, float>> allHits;
    df.Foreach([&](const RVec<short> &det, const RVec<short> &layer,
                   const RVec<float> &edge, const RVec<float> &chi2, const RVec<int16_t> &ndf,
                   const RVec<float> &theta, const RVec<int16_t> &pindex,
                   const RVec<int> &pid, const RVec<int> &passFid) {
        for (size_t i = 0; i < layer.size(); ++i) {
            int idx = pindex[i];
            if (det[i] == selecteddetector && pid[idx] == selectedPid && passFid[idx]) {
                float chi2ndf = ndf[idx] > 0 ? chi2[idx] / ndf[idx] : 0;
                allHits.emplace_back(layer[i], edge[i], chi2ndf, theta[idx]);
            }
        }
    }, {"REC_Traj_detector", "REC_Traj_layer", "REC_Traj_edge", "REC_Track_chi2", "REC_Track_NDF",
        "REC_Particle_theta", "REC_Traj_pindex", "REC_Particle_pid", "REC_Track_pass_fid"});

    std::map<std::string, TH2F*> histos;

    for (const auto &[layer, edge, chi2ndf, theta] : allHits) {
        auto it = std::find(layers.begin(), layers.end(), layer);
        if (it == layers.end()) continue;
        size_t li = std::distance(layers.begin(), it);
        int ti = thetaRegionIndex(theta, thetaCuts);
        std::string key = "L" + std::to_string(layer) + "_T" + std::to_string(ti);

        if (!histos.count(key)) {
            histos[key] = new TH2F(key.c_str(), ";edge [cm];<chi2/ndf>", nbinsx, xmins[li], xmaxs[li], 100, ymin, ymax);
            histos[key]->SetDirectory(nullptr);
        }
        histos[key]->Fill(edge, chi2ndf);
    }

    for (size_t li = 0; li < layers.size(); ++li) {
        TCanvas *c = new TCanvas(("c_layer" + std::to_string(layers[li])).c_str(), "", 2000, 1500);
        TLegend *leg = new TLegend(0.65, 0.7, 0.88, 0.88);
        leg->SetTextFont(42);
        leg->SetTextSize(0.03);
        //leg->SetHeader("#scale[1.0]{#bf{#theta ranges}}", "C");
        int color = 1;

        for (int ti = 0; ti <= (int)thetaCuts.size(); ++ti) {
            std::string key = "L" + std::to_string(layers[li]) + "_T" + std::to_string(ti);
            if (!histos.count(key)) continue;

            TProfile *pfx = histos[key]->ProfileX(("prof_" + key).c_str(), 1, -1, "s");
            pfx->SetBins(nbinsx, xmins[li], xmaxs[li]);
            pfx->GetYaxis()->SetRangeUser(0, 100);
            pfx->SetLineColor(color);
            pfx->SetMarkerColor(color);
            pfx->SetMarkerStyle(20);
            pfx->SetMarkerSize(1.0);
            pfx->SetLineWidth(2);
            pfx->SetTitle(Form("proton CVT layer %d", layers[li]));
            pfx->Draw(ti == 0 ? "PE" : "PE SAME");
            color++;

            std::string label;
            if (ti == 0)
                label = Form("#scale[0.9]{#theta/deg #leq %.0f}", thetaCuts[0]);
            else if (ti == (int)thetaCuts.size())
                label = Form("#scale[0.9]{#theta/deg > %.0f}", thetaCuts.back());
            else
                label = Form("#scale[0.9]{%.0f < #theta/deg #leq %.0f}", thetaCuts[ti - 1], thetaCuts[ti]);

            leg->AddEntry(pfx, label.c_str(), "lp");
        }
        leg->Draw();

        std::string outname = "proton_CVT_layer" + std::to_string(layers[li]) + ".png";
        c->SaveAs(outname.c_str());
        std::cout << "Saved: " << outname << std::endl;
        delete c;
    }
}

void analysisCVTFid() {
    std::string path = "./../build/";
    std::vector<int> layers = {1, 3, 5, 7, 12};
    std::vector<float> xmins = {-0.5, -0.5, -0.5, -4.0, -5.0};
    std::vector<float> xmaxs = {2.5, 2.5, 2.5, 20.0, 25.0};
    std::vector<float> thetaCuts = {55, 85, 115};
    DrawCVTChi2ndf_Optimized(2212, 5, 0, 400, 50, layers, xmins, xmaxs, thetaCuts, path + "dfSelected_after_fiducialCuts.root", "dfSelected_after");
    gApplication->Terminate(0);
}
