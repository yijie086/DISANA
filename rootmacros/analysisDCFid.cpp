#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TCanvas.h>
#include <TFile.h>
#include <TLegend.h>
#include <TMath.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TStopwatch.h>
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

void DrawDCChi2ndf_Optimized(const int &selectedPid, const int &selecteddetector,
                               const double ymin, const double ymax, const int nbinsx,
                               const std::vector<int> &layers,
                               const std::vector<int> &sectors,
                               const std::vector<float> &xmins,
                               const std::vector<float> &xmaxs,
                               const std::vector<float> &thetaCuts,
                               const std::string &filename, const std::string &treename, const bool doFid) {
    TStopwatch timer;
    timer.Start();
    ROOT::EnableImplicitMT();  // Enable multi-threading for RDataFrame
    ROOT::RDataFrame df(treename, filename);
    gStyle->SetOptStat(0);

    std::map<std::string, TH2F*> histos;

    for (size_t li = 0; li < layers.size(); ++li) {
        for (size_t ti = 0; ti <= thetaCuts.size(); ++ti) {
            for (size_t si = 0; si < sectors.size(); ++si) {
                std::string key = "L" + std::to_string(layers[li]) + "_S" + std::to_string(sectors[si]) + "_T" + std::to_string(ti);
                histos[key] = new TH2F(key.c_str(), ";edge [cm];<chi2/ndf>", nbinsx, xmins[li], xmaxs[li], 500, ymin, ymax);
                histos[key]->SetDirectory(nullptr);
            }
        }
    }

    df.Foreach([&](const RVec<short> &det, const RVec<short> &layer,
                    const RVec<float> &edge, const RVec<float> &theta,
                    const RVec<float> &chi2, const RVec<int16_t> &ndf,
                    const RVec<int16_t> &trackpindex, const RVec<short> &sector,
                    const RVec<int16_t> &pindex, const RVec<int> &pid, const RVec<int> &passFid) {
        for (size_t i = 0; i < layer.size(); ++i) {
            int idx = pindex[i];
            int fidpass = doFid ? passFid[idx] : 1;  // 如果不需要fiducial cuts，则始终通过
            if (det[i] != selecteddetector || pid[idx] != selectedPid || !fidpass) continue;
            float trackchi2 = 0;
            float trackndf = 0;
            short tracksector = 0;
            for (size_t j = 0; j < trackpindex.size(); ++j) {
                if (trackpindex[j] == idx) {
                    trackchi2 = chi2[j];
                    trackndf = ndf[j];
                    tracksector = sector[j];
                    break;  // 找到对应的trackpindex后退出循环
                }
            }
            float chi2ndf = trackndf > 0 ? trackchi2 / trackndf : 0;
            int lay = layer[i];
            float edg = edge[i];
            float th = theta[i];

            auto it = std::find(layers.begin(), layers.end(), lay);
            if (it == layers.end()) continue;
            size_t li = std::distance(layers.begin(), it);
            int ti = thetaRegionIndex(th, thetaCuts);
            if (tracksector < 1 || tracksector > 6) {
                std::cout << "Warning: tracksector is wrong for pindex " << idx << ", skipping this entry." << std::endl;
            }

            std::string key = "L" + std::to_string(lay) + "_S" + std::to_string(tracksector) + "_T" + std::to_string(ti);
            if (!histos.count(key)) {
                std::cerr << "Warning: histogram key not found -> " << key << std::endl;
                std::cout << key << " does not exist in histos map." << std::endl;
            }

            histos[key]->Fill(edg, chi2ndf);
        }
    }, {"REC_Traj_detector", "REC_Traj_layer", "REC_Traj_edge", "REC_Particle_theta", "REC_Track_chi2",
        "REC_Track_NDF", "REC_Track_pindex", "REC_Track_sector", "REC_Traj_pindex", "REC_Particle_pid", "REC_Track_pass_fid"});

    for (size_t si = 0; si < sectors.size(); ++si) {
        for (size_t li = 0; li < layers.size(); ++li) {
            TCanvas *c = new TCanvas(("c_layer" + std::to_string(layers[li]) + "_sector" + std::to_string(sectors[si])).c_str(), "", 2000, 1500);
            TLegend *leg = new TLegend(0.65, 0.7, 0.88, 0.88);
            leg->SetTextFont(42);
            leg->SetTextSize(0.03);
            std::vector<int> colorList = {
                kBlack,      // 1
                kRed+1,      // 2
                kBlue+1,     // 3
                kGreen+2,    // 4
                kMagenta+2,  // 5
                kCyan+1,     // 6
                kOrange+7,   // 7
                kViolet+1,   // 8
                kTeal+1,     // 9
                kAzure+1,    // 10
                kPink+1,     // 11
                kSpring+4,   // 12
                kYellow+2,   // 13
                kGray+2,     // 14
                kRed-7       // 15: 暗红
            };


            for (int ti = 0; ti <= (int)thetaCuts.size(); ++ti) {
                std::string key = "L" + std::to_string(layers[li]) + "_S" + std::to_string(sectors[si]) + "_T" + std::to_string(ti);
                if (!histos.count(key)) continue;

                TProfile *pfx = histos[key]->ProfileX(("prof_" + key).c_str(), 1, -1, "s");
                pfx->SetBins(nbinsx, xmins[li], xmaxs[li]);
                pfx->GetYaxis()->SetRangeUser(ymin, ymax/2);
                int color = colorList[ti % colorList.size()];
                pfx->SetLineColor(color);
                pfx->SetMarkerColor(color);
                pfx->SetMarkerStyle(20);
                pfx->SetMarkerSize(1.0);
                pfx->SetLineWidth(2);
                if (selectedPid == 2212){
                    if (doFid){
                        pfx->SetTitle(Form("proton DC layer %d sector %d", layers[li], sectors[si]));
                    } else {
                        pfx->SetTitle(Form("proton DC layer %d sector %d (no fiducial cuts)", layers[li], sectors[si]));
                    }
                }
                else if (selectedPid == 11) {
                    if (doFid){
                        pfx->SetTitle(Form("electron DC layer %d sector %d", layers[li], sectors[si]));
                    } else {
                        pfx->SetTitle(Form("electron DC layer %d sector %d (no fiducial cuts)", layers[li], sectors[si]));
                    }
                }
                else {
                    if (doFid){
                        pfx->SetTitle(Form("pid%d DC layer %d sector %d", selectedPid, layers[li], sectors[si]));
                    } else {
                        pfx->SetTitle(Form("pid%d DC layer %d sector %d (no fiducial cuts)", selectedPid, layers[li], sectors[si]));
                    }
                }

                //pfx->SetTitle(Form("proton CVT layer %d", layers[li]));
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

            std::string outname;
            if (selectedPid == 2212) {
                if (doFid)
                    outname = Form("proton_DC_layer%d_sector%d.png", layers[li], sectors[si]);
                else
                    outname = Form("proton_DC_layer%d_sector%d_noFid.png", layers[li], sectors[si]);
            }
            else if (selectedPid == 11) {
                if (doFid)
                    outname = Form("electron_DC_layer%d_sector%d.png", layers[li], sectors[si]);
                else
                    outname = Form("electron_DC_layer%d_sector%d_noFid.png", layers[li], sectors[si]);
            }
            else {
                if (doFid)
                    outname = Form("pid%d_DC_layer%d_sector%d.png", selectedPid, layers[li], sectors[si]);
                else
                    outname = Form("pid%d_DC_layer%d_sector%d_noFid.png", selectedPid, layers[li], sectors[si]);
            }
            std::string outdir = "DCChi2ndfPlots";
            gSystem->Exec(("mkdir -p " + outdir).c_str());
            outname = outdir + "/" + outname;
            c->SaveAs(outname.c_str());
            std::cout << "Saved: " << outname << std::endl;
            delete c;
        }
    }
    timer.Stop();
    std::cout << "Time for DrawDCChi2ndf_Optimized: ";
    timer.Print();
}

void DrawDCHitResponse(const int &selectedPid, const int &selecteddetector,
                        const std::vector<int> &layers,
                        const std::string &filename, const std::string &treename,
                        const bool doFid) {
    TStopwatch timer;
    timer.Start();
    ROOT::EnableImplicitMT();  // Enable multi-threading for RDataFrame
    ROOT::RDataFrame df(treename, filename);
    gStyle->SetOptStat(0);

    // 初始化所有层的TH2F
    std::map<int, TH2F*> histos;
    for (int layer : layers) {
        std::string name = Form("theta_vs_phi_L%d", layer);
        std::string title;
        if (selectedPid==2212){
            if (doFid) {
                name = Form("proton_theta_vs_phi_L%d", layer);
                title = Form("proton DC layer %d;#phi [deg];#theta [deg]", layer);
            } else {
                name = Form("proton_theta_vs_phi_L%d_noFid", layer);
                title = Form("proton DC layer %d (no fiducial cuts);#phi [deg];#theta [deg]", layer);
            }
        } else if (selectedPid == 11) {
            if (doFid) {
                name = Form("electron_theta_vs_phi_L%d", layer);
                title = Form("electron DC layer %d;#phi [deg];#theta [deg]", layer);
            } else {
                name = Form("electron_theta_vs_phi_L%d_noFid", layer);
                title = Form("electron DC layer %d (no fiducial cuts);#phi [deg];#theta [deg]", layer);
            }
        } else {
            if (doFid) {
                name = Form("pid%d_theta_vs_phi_L%d", selectedPid, layer);
                title = Form("pid%d DC layer %d;#phi [deg];#theta [deg]", selectedPid, layer);
            } else {
                name = Form("pid%d_theta_vs_phi_L%d_noFid", selectedPid, layer);
                title = Form("pid%d DC layer %d (no fiducial cuts);#phi [deg];#theta [deg]", selectedPid, layer);
            }
        }
        
        histos[layer] = new TH2F(name.c_str(), title.c_str(), 500, -350, 350, 500, -350, 350);
        histos[layer]->SetDirectory(nullptr);
    }

    df.Foreach([&](const RVec<short> &det, const RVec<short> &layer,
                             const RVec<float> &theta_deg, const RVec<float> &phi_deg,
                             const RVec<int16_t> &pindex, const RVec<int> &pid,
                             const RVec<int> &passFid) {
        for (size_t i = 0; i < det.size(); ++i) {
            if (det[i] != selecteddetector) continue;
            int lay = layer[i];
            int idx = pindex[i];
            int fidpass = doFid ? passFid[idx] : 1;  // 如果不需要fiducial cuts，则始终通过
            if (pid[idx] != selectedPid || !fidpass) continue;
            if (histos.count(lay))
                histos[lay]->Fill(phi_deg[i], theta_deg[i]);
        }
    }, {"REC_Traj_detector", "REC_Traj_layer", "REC_Traj_x", "REC_Traj_y",  "REC_Traj_pindex", "REC_Particle_pid", "REC_Track_pass_fid"});

    std::string outdir = "DCHitResponsePlots";
    gSystem->Exec(("mkdir -p " + outdir).c_str());

    for (int layer : layers) {
        TCanvas *c = new TCanvas(Form("c_hit_layer_%d", layer), "", 2700, 2000);
        gPad->SetLogz();
        histos[layer]->Draw("COLZ");

        std::string outname;
        if (selectedPid == 2212) {
            if (doFid)
                outname = Form("proton_DC_hit_layer%d.png", layer);
            else
                outname = Form("proton_DC_hit_layer%d_noFid.png", layer);
        } else if (selectedPid == 11) {
            if (doFid)
                outname = Form("electron_DC_hit_layer%d.png", layer);
            else
                outname = Form("electron_DC_hit_layer%d_noFid.png", layer);
        } else {
            if (doFid)
                outname = Form("pid%d_DC_hit_layer%d.png", selectedPid, layer);
            else
                outname = Form("pid%d_DC_hit_layer%d_noFid.png", selectedPid, layer);
        }
        outname = outdir + "/" + outname;

        c->SaveAs(outname.c_str());
        std::cout << "Saved: " << outname << std::endl;
        delete c;
    }


    timer.Stop();
    std::cout << "Time for DrawCVTHitResponse: ";
    timer.Print();
}


void analysisDCFid() {
    //std::string path = "/work/clas12/yijie/clas12ana/analysis203/DISANA/build/bbbs/";
    std::string path = "./../build/";
    std::vector<int> layers = {6, 18, 36};
    std::vector<float> xmins = {0, 0, 0, 0, 0, 0};
    std::vector<float> xmaxs = {25, 25, 25, 25, 25, 25};
    std::vector<float> thetaCuts = {10, 15, 20, 25};
    std::vector<int> sectors = {1, 2, 3, 4, 5, 6};
    DrawDCChi2ndf_Optimized(11, 6, 0, 400, 50, layers, sectors, xmins, xmaxs, thetaCuts, path + "dfSelected_after_fiducialCuts.root", "dfSelected_after", true);
    DrawDCChi2ndf_Optimized(11, 6, 0, 400, 50, layers, sectors, xmins, xmaxs, thetaCuts, path + "dfSelected_before_fiducialCuts.root", "dfSelected_before", false);
    DrawDCChi2ndf_Optimized(2212, 6, 0, 400, 50, layers, sectors, xmins, xmaxs, thetaCuts, path + "dfSelected_after_fiducialCuts.root", "dfSelected_after", true);
    DrawDCChi2ndf_Optimized(2212, 6, 0, 400, 50, layers, sectors, xmins, xmaxs, thetaCuts, path + "dfSelected_before_fiducialCuts.root", "dfSelected_before", false);
    DrawDCHitResponse(11, 6, layers, path + "dfSelected_after_fiducialCuts.root", "dfSelected_after",true);
    DrawDCHitResponse(2212, 6, layers, path + "dfSelected_after_fiducialCuts.root", "dfSelected_after",true);
    DrawDCHitResponse(11, 6, layers, path + "dfSelected_before_fiducialCuts.root", "dfSelected_before",false);
    DrawDCHitResponse(2212, 6, layers, path + "dfSelected_before_fiducialCuts.root", "dfSelected_before",false);
    gApplication->Terminate(0);
}
