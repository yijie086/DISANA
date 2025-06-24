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

void DrawECALHitResponse(const int &selectedPid, const int &selecteddetector,
                         const std::vector<int> &layers,
                         const std::vector<int> &sectors,
                         const std::string &filename, const std::string &treename,
                         const bool doFid) {
    TStopwatch timer;
    timer.Start();
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df(treename, filename);
    gStyle->SetOptStat(0);

    std::map<std::string, TH2F*> histos;

    for (int layer : layers) {
        for (int sector : sectors) {
            std::string name, title;
            std::string key = Form("pid%d_ECAL_L%d_S%d%s", selectedPid, layer, sector, doFid ? "" : "_noFid");
            if (selectedPid == 2212) {
                if (doFid) {
                    name = Form("proton_ECAL_L%d_S%d", layer, sector);
                    title = Form("proton ECAL layer %d sector %d", layer, sector);
                } else {
                    name = Form("proton_ECAL_L%d_S%d_noFid", layer, sector);
                    title = Form("proton ECAL layer %d sector %d (no fiducial cuts)", layer, sector);
                }
            } else if (selectedPid == 11) {
                if (doFid) {
                    name = Form("electron_ECAL_L%d_S%d", layer, sector);
                    title = Form("electron ECAL layer %d sector %d", layer, sector);
                } else {
                    name = Form("electron_ECAL_L%d_S%d_noFid", layer, sector);
                    title = Form("electron ECAL layer %d sector %d (no fiducial cuts)", layer, sector);
                }
            } else {
                if (doFid) {
                    name = Form("pid%d_ECAL_L%d_S%d", selectedPid, layer, sector);
                    title = Form("pid%d ECAL layer %d sector %d", selectedPid, layer, sector);
                } else {
                    name = Form("pid%d_ECAL_L%d_S%d_noFid", selectedPid, layer, sector);
                    title = Form("pid%d ECAL layer %d sector %d (no fiducial cuts)", selectedPid, layer, sector);
                }
            }

            //title += ";x [cm];y [cm]";
            histos[key] = new TH2F(name.c_str(), title.c_str(), 500, 0, 500, 500, 0, 500);
            histos[key]->SetDirectory(nullptr);
        }
    }

    df.Foreach([&](const RVec<short> &det, const RVec<short> &layer,
                   const RVec<short> &sector, const RVec<int16_t> &pindex,
                   const RVec<int> &pid, const RVec<int> &passFid,
                   const RVec<float> &lu, const RVec<float> &lv, const RVec<float> &lw) {
        for (size_t i = 0; i < det.size(); ++i) {
            if (det[i] != selecteddetector) continue;
            int lay = layer[i];
            int sec = sector[i];
            int idx = pindex[i];
            int fidpass = doFid ? passFid[idx] : 1;
            if (pid[idx] != selectedPid || !fidpass) continue;
            std::string key = Form("pid%d_ECAL_L%d_S%d%s", selectedPid, lay, sec, doFid ? "" : "_noFid");
            if (!histos.count(key)) continue;

            float x = 0, y = 0;
            if (lay == 1 || lay == 4) {
                x = lw[i];
                y = lv[i];
                histos[key]->GetXaxis()->SetTitle("lw [cm]");
                histos[key]->GetYaxis()->SetTitle("lv [cm]");
            } else if (lay == 7) {
                x = lu[i];
                y = lv[i];
                histos[key]->GetXaxis()->SetTitle("lu [cm]");
                histos[key]->GetYaxis()->SetTitle("lv [cm]");
            }
            histos[key]->Fill(x, y);
        }
    }, {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector",
        "REC_Calorimeter_pindex", "REC_Particle_pid", "REC_Track_pass_fid",
        "REC_Calorimeter_lu", "REC_Calorimeter_lv", "REC_Calorimeter_lw"});

    std::string outdir = "ECALHitResponsePlots";
    gSystem->Exec(("mkdir -p " + outdir).c_str());

    for (int layer : layers) {
        for (int sector : sectors) {
            std::string key = Form("pid%d_ECAL_L%d_S%d%s", selectedPid, layer, sector, doFid ? "" : "_noFid");
            if (!histos.count(key)) continue;

            TCanvas *c = new TCanvas(Form("c_ecal_hit_layer%d_sector%d", layer, sector), "", 2700, 2000);
            gPad->SetLogz();
            histos[key]->Draw("COLZ");

            std::string outname = Form("%s/pid%d_ECAL_hit_layer%d_sector%d%s.png",
                                       outdir.c_str(), selectedPid, layer, sector, doFid ? "" : "_noFid");
            if (selectedPid == 11) {
                outname = Form("%s/electron_ECAL_hit_layer%d_sector%d%s.png",
                               outdir.c_str(), layer, sector, doFid ? "" : "_noFid");
            } else if (selectedPid == 22) {
                outname = Form("%s/photon%d_sector%d%s.png",
                               outdir.c_str(), layer, sector, doFid ? "" : "_noFid");
            } else {
                outname = Form("%s/pid%d_ECAL_hit_layer%d_sector%d%s.png",
                               outdir.c_str(), selectedPid, layer, sector, doFid ? "" : "_noFid");
            }
            c->SaveAs(outname.c_str());
            std::cout << "Saved: " << outname << std::endl;
            delete c;
        }
    }

    timer.Stop();
    std::cout << "Time for DrawECALHitResponse: ";
    timer.Print();
}


void DrawECALEnergyProfile(const int &selectedPid, const int &selecteddetector,
                                             const std::vector<int> &layers,
                                             const std::vector<int> &sectors,
                                             const std::string &filename, const std::string &treename,
                                             const bool doFid) {
    TStopwatch timer;
    timer.Start();
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df(treename, filename);
    gStyle->SetOptStat(0);

    std::map<int, TH2F*> hist_lw, hist_lv;

    for (int sector : sectors) {
        hist_lw[sector] = new TH2F(Form("lw_S%d", sector), ";PCAL lw [cm];E/p", 500, 0, 100, 500, 0, 2);
        hist_lv[sector] = new TH2F(Form("lv_S%d", sector), ";PCAL lv [cm];E/p", 500, 0, 100, 500, 0, 2);
        hist_lw[sector]->SetDirectory(nullptr);
        hist_lv[sector]->SetDirectory(nullptr);
    }

    df.Foreach([&](const RVec<short> &det, const RVec<short> &layer,
                   const RVec<short> &sector, const RVec<int16_t> &pindex,
                   const RVec<int> &pid, const RVec<int> &passFid,
                   const RVec<float> &lv, const RVec<float> &lw,
                   const RVec<float> &energy,
                   const RVec<float> &px, const RVec<float> &py, const RVec<float> &pz) {
        struct HitInfo {
            float totalE = 0;
            float lw_layer1 = -999;
            float lv_layer1 = -999;
            bool hasL1 = false;
        };

        std::map<std::pair<int16_t, int>, HitInfo> hitMap;  // key: (pindex, sector)

        for (size_t i = 0; i < det.size(); ++i) {
            if (det[i] != selecteddetector) continue;
            int lay = layer[i];
            int sec = sector[i];
            int16_t idx = pindex[i];
            if (pid[idx] != selectedPid || (doFid && !passFid[idx])) continue;
            if (std::find(layers.begin(), layers.end(), lay) == layers.end()) continue;

            auto key = std::make_pair(idx, sec);
            hitMap[key].totalE += energy[i];
            if (lay == 1) {
                hitMap[key].lw_layer1 = lw[i];
                hitMap[key].lv_layer1 = lv[i];
                hitMap[key].hasL1 = true;
            }
        }

        for (const auto &[key, info] : hitMap) {
            int16_t idx = key.first;
            int sec = key.second;
            if (!info.hasL1) continue;  // skip if no layer 1
            if (!hist_lw.count(sec)) continue;
            if (idx >= (int)px.size()) continue;

            float p = std::sqrt(px[idx]*px[idx] + py[idx]*py[idx] + pz[idx]*pz[idx]);
            if (p < 1e-3) continue;

            float eOverP = info.totalE / p;
            hist_lw[sec]->Fill(info.lw_layer1, eOverP);
            hist_lv[sec]->Fill(info.lv_layer1, eOverP);
        }
    }, {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector",
        "REC_Calorimeter_pindex", "REC_Particle_pid", "REC_Track_pass_fid",
        "REC_Calorimeter_lv", "REC_Calorimeter_lw", "REC_Calorimeter_energy",
        "REC_Particle_px", "REC_Particle_py", "REC_Particle_pz"});

    std::string outdir = "ECALdepositedEnergyPlots";
    gSystem->Exec(("mkdir -p " + outdir).c_str());

    for (int sector : sectors) {
        if (!hist_lw.count(sector)) continue;

        TCanvas *c = new TCanvas(Form("c_ecal_EoverP_S%d", sector), "", 2000, 1500);
        TLegend *leg = new TLegend(0.6, 0.75, 0.88, 0.88);
        leg->SetTextSize(0.035);
        leg->SetBorderSize(0);

        TProfile *p_lw = hist_lw[sector]->ProfileX(Form("pfx_lw_S%d", sector), 1, -1, "s");
        p_lw->SetLineColor(kBlue + 1);
        p_lw->SetMarkerColor(kBlue + 1);
        p_lw->SetMarkerStyle(20);
        p_lw->SetTitle(Form("E/p in Sector %d (Position from PCAL); lv or lw [cm];E/p", sector));
        p_lw->GetYaxis()->SetRangeUser(0.15, 0.3);
        p_lw->Draw("PE");
        leg->AddEntry(p_lw, "E/p vs lw", "lp");

        TProfile *p_lv = hist_lv[sector]->ProfileX(Form("pfx_lv_S%d", sector), 1, -1, "s");
        p_lv->SetLineColor(kRed + 1);
        p_lv->SetMarkerColor(kRed + 1);
        p_lv->SetMarkerStyle(21);
        p_lv->Draw("PE SAME");
        leg->AddEntry(p_lv, "E/p vs lv", "lp");

        leg->Draw();

        std::string outname;
        if (selectedPid == 11) {
            outname = Form("%s/electron_EoverP_L1_sector%d%s.png", outdir.c_str(), sector, doFid ? "" : "_noFid");
        } else if (selectedPid == 22) {
            outname = Form("%s/photon_EoverP_L1_sector%d%s.png", outdir.c_str(), sector, doFid ? "" : "_noFid");
        } else {
            outname = Form("%s/pid%d_EoverP_L1_sector%d%s.png", outdir.c_str(), selectedPid, sector, doFid ? "" : "_noFid");
        }

        c->SaveAs(outname.c_str());
        std::cout << "Saved: " << outname << std::endl;
        delete c;
    }

    timer.Stop();
    std::cout << "Time for DrawECALEnergyProfile: ";
    timer.Print();
}




void analysisECALFid() {
    //std::string path = "/work/clas12/yijie/clas12ana/analysis203/DISANA/build/bbbs/";
    //std::string path = "./../build/";
    std::string path = "/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/NewAnalysisFrameWork/testing_outupt/afterFiducialCuts/afterallFidCuts_dsts/";
    std::vector<int> layers = {1, 4, 7};
    std::vector<int> sectors = {1, 2, 3, 4, 5, 6};
    DrawECALHitResponse(11, 7, layers, sectors ,path + "dfSelected.root", "dfSelected",false);
    DrawECALHitResponse(11, 7, layers, sectors ,path + "dfSelected_afterFid.root", "dfSelected_afterFid",true);
    DrawECALEnergyProfile(11, 7, layers, sectors, path + "dfSelected_afterFid.root", "dfSelected_afterFid", true);
    DrawECALEnergyProfile(11, 7, layers, sectors, path + "dfSelected.root", "dfSelected", false);
    gApplication->Terminate(0);
}
