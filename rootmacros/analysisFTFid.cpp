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
                        const std::string &filename, const std::string &treename) {
    TStopwatch timer;
    timer.Start();

    ROOT::EnableImplicitMT();  // Enable multi-threading for RDataFrame

    ROOT::RDataFrame df(treename, filename);
    gStyle->SetOptStat(0);

    // 初始化所有层的TH2F
    std::map<int, TH2F*> histos;
    for (int layer : layers) {
        std::string name = Form("x_vs_y_L%d", layer);
        std::string title = Form("photon FTCal layer %d;x [cm];y [cm]", layer);
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
            if (pid[idx] != selectedPid || !passFid[idx]) continue;
            if (histos.count(lay))
                histos[lay]->Fill(x[i], y[i]);
        }
    }, {"REC_ForwardTagger_detector", "REC_ForwardTagger_layer", "REC_ForwardTagger_x", "REC_ForwardTagger_y",  "REC_ForwardTagger_pindex", "REC_Particle_pid", "REC_Track_pass_fid"});

    for (int layer : layers) {
        TCanvas *c = new TCanvas(Form("c_hit_layer_%d", layer), "", 2400, 2000);
        histos[layer]->Draw("COLZ");
        //gPad->SetLogz();

        c->SaveAs(Form("photon_FTCal_hit_layer%d.png", layer));
        delete c;
    }

    timer.Stop();
    std::cout << "Time for DrawCVTHitResponse: ";
    timer.Print();
}


void analysisFTFid() {
    //std::string path = "/work/clas12/yijie/clas12ana/analysis203/DISANA/build/bbbs/";
    std::string path = "./../build/";
    std::vector<int> layers = {1};
    DrawFTHitResponse(22, 10, layers, path + "dfSelected_before_fiducialCuts.root", "dfSelected_before");
    gApplication->Terminate(0);
}
