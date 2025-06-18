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

void DrawFTHitResponse(const int &selectedPid, const int &selecteddetector,
                        const std::vector<int> &layers,
                        const std::string &filename, const std::string &treename, const bool doFid) {
    TStopwatch timer;
    timer.Start();

    ROOT::EnableImplicitMT();  // Enable multi-threading for RDataFrame

    ROOT::RDataFrame df(treename, filename);
    gStyle->SetOptStat(0);

    // 初始化所有层的TH2F
    std::map<int, TH2F*> histos;
    for (int layer : layers) {
        std::string name;
        std::string title;
        if (selectedPid == 22) {  // Photon
            if (doFid) {
                name = Form("photon_x_vs_y_L%d", layer);
                title = Form("photon FTCal layer %d;x [cm];y [cm]", layer);
            } else {
                name = Form("photon_x_vs_y_L%d_noFid", layer);
                title = Form("photon FTCal layer %d (no fiducial cuts);x [cm];y [cm]", layer);
            }
        }
        else if (selectedPid == 11) {  // Electron
            if (doFid) {
                name = Form("electron_x_vs_y_L%d", layer);
                title = Form("electron FTCal layer %d;x [cm];y [cm]", layer);
            } else {
                name = Form("electron_x_vs_y_L%d_noFid", layer);
                title = Form("electron FTCal layer %d (no fiducial cuts);x [cm];y [cm]", layer);
            }
        }
        else {
            if (doFid) {
                name = Form("pid%d_x_vs_y_L%d", selectedPid, layer);
                title = Form("pid%d FTCal layer %d;x [cm];y [cm]", selectedPid, layer);
            } else {
                name = Form("pid%d_x_vs_y_L%d_noFid", selectedPid, layer);
                title = Form("pid%d FTCal layer %d (no fiducial cuts);x [cm];y [cm]", selectedPid, layer);
            }
        }
        histos[layer] = new TH2F(name.c_str(), title.c_str(), 400, -20, 20, 400, -20, 20);
        histos[layer]->SetDirectory(nullptr);
    }

    df.Foreach([&](const RVec<short> &det, const RVec<short> &layer,
                             const RVec<float> &x, const RVec<float> &y,
                             const RVec<int16_t> &pindex, const RVec<int> &pid,
                             const RVec<int> &passFid) {
        for (size_t i = 0; i < det.size(); ++i) {
            if (det[i] != selecteddetector) continue;
            int lay = layer[i];
            int idx = pindex[i];
            int fidpass = doFid ? passFid[idx] : 1;  // 如果不需要fiducial cuts，则始终通过
            if (pid[idx] != selectedPid || !fidpass) continue;
            if (histos.count(lay))
                histos[lay]->Fill(x[i], y[i]);
        }
    }, {"REC_ForwardTagger_detector", "REC_ForwardTagger_layer", "REC_ForwardTagger_x", "REC_ForwardTagger_y",  "REC_ForwardTagger_pindex", "REC_Particle_pid", "REC_Track_pass_fid"});

    std::string outdir = "FTCalHitPlots";
    gSystem->Exec(("mkdir -p " + outdir).c_str());

    for (int layer : layers) {
        TCanvas *c = new TCanvas(Form("c_hit_layer_%d", layer), "", 2400, 2000);
        histos[layer]->Draw("COLZ");

        std::string outname;
        if (selectedPid == 22) {
            outname = outdir + "/" + Form(doFid ? "photon_FTCal_hit_layer%d.png" : "photon_FTCal_hit_layer%d_noFid.png", layer);
        } else if (selectedPid == 11) {
            outname = outdir + "/" + Form(doFid ? "electron_FTCal_hit_layer%d.png" : "electron_FTCal_hit_layer%d_noFid.png", layer);
        } else {
            outname = outdir + "/" + Form(doFid ? "pid%d_FTCal_hit_layer%d.png" : "pid%d_FTCal_hit_layer%d_noFid.png", selectedPid, layer);
        }

        c->SaveAs(outname.c_str());
        std::cout << "Saved: " << outname << std::endl;
        delete c;
    }

    timer.Stop();
    std::cout << "Time for DrawFTHitResponse: ";
    timer.Print();
}


void analysisFTFid() {
    //std::string path = "/work/clas12/yijie/clas12ana/analysis203/DISANA/build/bbbs/";
    std::string path = "./../build/";
    std::vector<int> layers = {1};
    DrawFTHitResponse(22, 10, layers, path + "dfSelected_after_fiducialCuts.root", "dfSelected_after",true);
    DrawFTHitResponse(22, 10, layers, path + "dfSelected_before_fiducialCuts.root", "dfSelected_before",false);
    DrawFTHitResponse(11, 10, layers, path + "dfSelected_after_fiducialCuts.root", "dfSelected_after",true);
    DrawFTHitResponse(11, 10, layers, path + "dfSelected_before_fiducialCuts.root", "dfSelected_before",false);
    gApplication->Terminate(0);
}
