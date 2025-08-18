#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TStopwatch.h>
#include <TSystem.h>
#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <map>
#include <cmath>
#include <TF1.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TSpline.h>
#include <memory>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>

using namespace ROOT::VecOps;

//================ Utility =================
int GetThetaRegionIndex(float thetaDeg, const std::vector<float> &thetaCuts) {
    for (size_t i = 0; i < thetaCuts.size(); ++i)
        if (thetaDeg <= thetaCuts[i]) return i;
    return thetaCuts.size();
}

std::string GetThetaBinLabel(int ti, const std::vector<float>& thetaCuts) {
    char buf[64];
    if (ti == 0)
        sprintf(buf, "#theta #leq %.1f#circ", thetaCuts[0]);
    else if (ti == (int)thetaCuts.size())
        sprintf(buf, "#theta > %.1f#circ", thetaCuts.back());
    else
        sprintf(buf, "%.1f#circ < #theta #leq %.1f#circ", thetaCuts[ti-1], thetaCuts[ti]);
    return std::string(buf);
}

std::string GetDetectorPart(short status) {
    short absStatus = std::abs(status);
    if (absStatus >= 1000 && absStatus < 2000) return "FT";
    else if (absStatus >= 2000 && absStatus < 3000)  return "FD";
    else if (absStatus >= 4000 && absStatus < 5000) return "CD";
    else return "Unknown";
}

TGraph* MakePeakGraph(TH2D* hist, const std::string& outTxtPath)
{
    const int minEntries = 100;
    const double maxChi2Ndf = 100.0;
    const double maxerrPeak = 0.2;

    std::vector<double> xVec, yVec, yerrVec, sigmaVec, chi2ndfVec;

    int nxbins = hist->GetNbinsX();
    for (int ix = 1; ix <= nxbins; ++ix) {
        std::unique_ptr<TH1D> proj(hist->ProjectionY(Form("proj_y_%d", ix), ix, ix, "e"));
        if (!proj || proj->GetEntries() < minEntries) continue;

        int maxBin = proj->GetMaximumBin();
        double mu0 = proj->GetBinCenter(maxBin);
        double sigma0 = proj->GetRMS();
        if (sigma0 == 0) continue;

        TF1 gaus("gaus", "gaus", mu0 - 2*sigma0, mu0 + 2*sigma0);
        gaus.SetParameters(proj->GetMaximum(), mu0, sigma0);
        int fitStatus = proj->Fit(&gaus, "QNR");

        bool fitOK = (fitStatus == 0) &&
                     (gaus.GetParameter(2) > 0) &&
                     (std::abs(gaus.GetParError(1)) < maxerrPeak) &&
                     (gaus.GetChisquare() / std::max(1.0, static_cast<double>(gaus.GetNDF())) < maxChi2Ndf);

        if (!fitOK) continue;

        double xC = hist->GetXaxis()->GetBinCenter(ix);
        double yPeak = gaus.GetParameter(1);
        double peakError = gaus.GetParError(1);
        double sigma = gaus.GetParameter(2);
        double chi2ndf = gaus.GetChisquare() / std::max(1.0, static_cast<double>(gaus.GetNDF()));

        xVec.push_back(xC);
        yVec.push_back(yPeak);
        yerrVec.push_back(peakError);
        sigmaVec.push_back(sigma);
        chi2ndfVec.push_back(chi2ndf);
    }

    // 删除边界点
    if (xVec.size() > 3) {
        xVec.erase(xVec.begin());        xVec.erase(xVec.begin());
        yVec.erase(yVec.begin());        yVec.erase(yVec.begin());
        yerrVec.erase(yerrVec.begin()); yerrVec.erase(yerrVec.begin());
        sigmaVec.erase(sigmaVec.begin()); sigmaVec.erase(sigmaVec.begin());
        chi2ndfVec.erase(chi2ndfVec.begin()); chi2ndfVec.erase(chi2ndfVec.begin());

        xVec.pop_back();    yVec.pop_back();
        yerrVec.pop_back();
        sigmaVec.pop_back(); chi2ndfVec.pop_back();
    }


    // 写入 TXT 文件
    std::ofstream fout(outTxtPath);
    fout << "# xCenter  yPeak yerr sigma  chi2/ndf\n";
    for (size_t i = 0; i < xVec.size(); ++i)
        fout << xVec[i] << "\t" << yVec[i] << "\t" << yerrVec[i]<<"\t" << sigmaVec[i] << "\t" << chi2ndfVec[i] << "\n";
    fout.close();
/*
    // 生成 TGraph
    TGraph* g = new TGraph(xVec.size(), xVec.data(), yVec.data());
    g->SetLineColor(kRed);
    g->SetLineWidth(3);
    return g;
*/
    auto *gErr = new TGraphErrors(xVec.size());
    for (size_t i = 0; i < xVec.size(); ++i) {
        double ex = 0.5 * hist->GetXaxis()->GetBinWidth(               // Bin half-width
                         hist->GetXaxis()->FindBin(xVec[i]));
        double ey = yerrVec[i];                                       // σ from fit
        gErr->SetPoint(i, xVec[i], yVec[i]);           // (xi , yi)
        gErr->SetPointError(i, ex, ey);                // (Δx , Δy)
    }

    // 样式：点 + 误差条
    gErr->SetMarkerStyle(20);
    gErr->SetMarkerSize(0.9);
    gErr->SetLineWidth(1);
    gErr->SetLineColor(kRed);
    gErr->SetMarkerColor(kRed);

    return gErr;
}


//================ Core =================
void DrawParticleKinematicsByThetaBins(
    const int &selectedPid,
    const std::vector<float> &thetaCuts,
    const std::string selecteddetector,                                        // in degree
    const std::vector<std::tuple<std::string,std::string,std::string,int,double,double>> &plotVars,     // {name, nbins, xmin, xmax}
    const std::string &filename,
    const std::string &treename)
    {
    TStopwatch timer;
    timer.Start();

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df(treename, filename);
    //gStyle->SetOptStat(1110);

    struct VarInfo {
        std::string saveName;
        std::string title;  // for future use
        std::string name;
        int nbins;
        double xmin, xmax;
        TH1D* overall;
        std::vector<TH1D*> binHists;   // size = thetaCuts.size()+1
    };

    //==== create histograms ====
    std::vector<VarInfo> vars;
    for (auto &cfg : plotVars) {
        VarInfo v;
        v.saveName = std::get<0>(cfg);
        v.title    = std::get<1>(cfg);
        v.name  = std::get<2>(cfg);
        v.nbins = std::get<3>(cfg);
        v.xmin  = std::get<4>(cfg);
        v.xmax  = std::get<5>(cfg);
        v.overall = new TH1D((v.name+"_overall").c_str(), "", v.nbins, v.xmin, v.xmax);
        v.overall->SetDirectory(nullptr);
        v.binHists.resize(thetaCuts.size()+1);
        for (size_t ti=0; ti<=thetaCuts.size(); ++ti) {
            std::string hname = v.name + Form("_T%zu", ti);
            v.binHists[ti] = new TH1D(hname.c_str(), "", v.nbins, v.xmin, v.xmax);
            v.binHists[ti]->SetDirectory(nullptr);
        }
        vars.push_back(v);
    }

    //==== fill histograms ====
    std::vector<std::string> colNames = {"REC_Particle_pid","REC_Particle_theta","REC_Particle_phi","REC_Particle_p","REC_Particle_status","REC_Particle_pass","REC_Particle_vz","REC_Particle_beta"};

    df.Foreach([&](const RVec<int> &pid,
                   const RVec<float> &theta,
                   const RVec<float> &phi,
                   const RVec<float> &p,
                   const RVec<short> &status,
                   const RVec<bool> &passPar,
                   const RVec<float> &vz,
                   const RVec<float> &beta
                   )  // passFid is not used but included for consistency
    {
        for (size_t i=0;i<pid.size();++i) {
            if (pid[i]!=selectedPid) continue;
            if (!passPar[i]) continue;  // only fill if the particle passes the fiducial cuts
            if (GetDetectorPart(status[i])!=selecteddetector && selecteddetector!="ALL") continue;  // check detector part
            float thetaDeg = theta[i]*180.0/M_PI;
            int ti = GetThetaRegionIndex(thetaDeg, thetaCuts);
            for (auto &v: vars) {
                double value = NAN;
                if (v.name=="theta") value = theta[i]*180.0/M_PI;
                else if (v.name=="phi"){value = phi[i]*180.0/M_PI; if (value < 0) value += 360.0;} // ensure phi is in [0, 360)
                else if (v.name=="p") value = p[i];
                else if (v.name=="vz") value = vz[i];  // assuming vz is the same as p for this example
                else if (v.name=="beta") value = beta[i];
                if (std::isnan(value)) continue;
                v.binHists[ti]->Fill(value);
                v.overall->Fill(value);
            }
        }
    }, colNames);

    //==== output directory ====
    std::string outDir = "ParticleKinematicPlots";
    gSystem->Exec(("mkdir -p "+outDir).c_str());

    //==== draw and save ====
    for (auto &v: vars) {
        // overall
        {
            TCanvas *c = new TCanvas(("c_"+v.name+"_overall").c_str(),"",1600,1200);
            v.overall->SetTitle((selecteddetector + " " +v.title + " " + v.name).c_str());
            if (v.name == "phi" || v.name == "theta") v.overall->GetXaxis()->SetTitle((v.name+" [deg]").c_str());
            else v.overall->GetXaxis()->SetTitle((v.name+" GeV").c_str());
            v.overall->GetYaxis()->SetTitle("Counts");
            v.overall->Draw();
            std::string out = outDir + "/" + v.saveName + "_" +selecteddetector+ "_" + v.name + "_overall.png";
            c->SaveAs(out.c_str());
            delete c;
            std::cout << "Saved: " << out << std::endl;
        }
        // theta bins
        for (size_t ti=0; ti<=thetaCuts.size(); ++ti) {
            TCanvas *c = new TCanvas(("c_"+v.name+Form("_T%zu",ti)).c_str(),"",1600,1200);
            if (v.name == "phi" || v.name == "theta") v.binHists[ti]->GetXaxis()->SetTitle((v.name+ " [deg]").c_str());
            else v.binHists[ti]->GetXaxis()->SetTitle((v.name+" [GeV]").c_str()); 
            v.binHists[ti]->SetTitle((selecteddetector + " " + v.title+" "+v.name + " in " + GetThetaBinLabel(ti, thetaCuts)).c_str());
            v.binHists[ti]->GetYaxis()->SetTitle("Counts");
            v.binHists[ti]->Draw();
            std::string out = outDir + "/" + v.saveName + "_" + selecteddetector + "_" + v.name + Form("_T%zu.png",ti);
            c->SaveAs(out.c_str());
            delete c;
            std::cout << "Saved: " << out << std::endl;
        }
    }

    timer.Stop();
    std::cout << "Time for DrawParticleKinematicsByThetaBins: ";
    timer.Print();
}

void Draw2DParticleKinematicsByThetaBins(
    const int &selectedPid,
    const std::vector<float> &thetaCuts,
    const std::string selecteddetector,
    const std::vector<std::tuple<std::string,std::string,std::string,int,double,double,int,double,double>> &plotVars,
    const std::string &filename,
    const std::string &treename) {
    TStopwatch timer;
    timer.Start();

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df(treename, filename);

    struct Var2DInfo {
        std::string saveName;
        std::string title;
        std::string name; // format: x:y (e.g., "p:theta")
        int nxbins;
        double xmin, xmax;
        int nybins;
        double ymin, ymax;
        TH2D* overall;
        std::vector<TH2D*> binHists;
    };

    std::vector<Var2DInfo> vars;
    for (auto &cfg : plotVars) {
        Var2DInfo v;
        v.saveName = std::get<0>(cfg);
        v.title    = std::get<1>(cfg);
        v.name     = std::get<2>(cfg);
        v.nxbins   = std::get<3>(cfg);
        v.xmin     = std::get<4>(cfg);
        v.xmax     = std::get<5>(cfg);
        v.nybins   = std::get<6>(cfg);
        v.ymin     = std::get<7>(cfg);
        v.ymax     = std::get<8>(cfg);

        v.overall = new TH2D((v.saveName+"_overall").c_str(), "", v.nxbins, v.xmin, v.xmax, v.nybins, v.ymin, v.ymax);
        v.overall->SetDirectory(nullptr);

        v.binHists.resize(thetaCuts.size() + 1);
        for (size_t ti = 0; ti <= thetaCuts.size(); ++ti) {
            std::string hname = v.saveName + Form("_T%zu", ti);
            v.binHists[ti] = new TH2D(hname.c_str(), "", v.nxbins, v.xmin, v.xmax, v.nybins, v.ymin, v.ymax);
            v.binHists[ti]->SetDirectory(nullptr);
        }
        vars.push_back(v);
    }

    std::vector<std::string> colNames = {"REC_Particle_pid", "REC_Particle_theta", "REC_Particle_phi", "REC_Particle_p", "REC_Particle_status", "REC_Particle_pass"};

    df.Foreach([&](const RVec<int> &pid,
                   const RVec<float> &theta,
                   const RVec<float> &phi,
                   const RVec<float> &p,
                   const RVec<short> &status,
                   const RVec<bool> &passPar)
    {
        for (size_t i = 0; i < pid.size(); ++i) {
            if (pid[i] != selectedPid) continue;
            if (!passPar[i]) continue;
            if (GetDetectorPart(status[i]) != selecteddetector && selecteddetector != "ALL") continue;

            float thetaDeg = theta[i] * 180.0 / M_PI;
            float phiDeg = phi[i] * 180.0 / M_PI;
            if (phiDeg < 0) phiDeg += 360.0;
            int ti = GetThetaRegionIndex(thetaDeg, thetaCuts);

            for (auto &v : vars) {
                double xval = NAN, yval = NAN;
                if (v.name == "p:theta") { xval = p[i]; yval = thetaDeg; }
                else if (v.name == "theta:p") { xval = thetaDeg; yval = p[i]; }
                else if (v.name == "phi:theta") { xval = phiDeg; yval = thetaDeg; }
                else if (v.name == "theta:phi") { xval = thetaDeg; yval = phiDeg; }
                else continue;

                if (std::isnan(xval) || std::isnan(yval)) continue;
                v.binHists[ti]->Fill(xval, yval);
                v.overall->Fill(xval, yval);
            }
        }
    }, colNames);

    std::string outDir = "Particle2DKinematicPlots";
    gSystem->Exec(("mkdir -p " + outDir).c_str());

    for (auto &v : vars) {
        {
            TCanvas *c = new TCanvas(("c_" + v.saveName + "_overall").c_str(), "", 1600, 1200);
            v.overall->SetTitle((selecteddetector + " " + v.title).c_str());
            v.overall->Draw("COLZ");
            if (v.name == "p:theta") {
                v.overall->GetXaxis()->SetTitle("p [GeV]");
                v.overall->GetYaxis()->SetTitle("#theta [deg]");
            } else if (v.name == "theta:p") {
                v.overall->GetXaxis()->SetTitle("#theta [deg]");
                v.overall->GetYaxis()->SetTitle("p [GeV]");
            } else if (v.name == "phi:theta") {
                v.overall->GetXaxis()->SetTitle("#phi [deg]");
                v.overall->GetYaxis()->SetTitle("#theta [deg]");
            } else if (v.name == "theta:phi") {
                v.overall->GetXaxis()->SetTitle("#theta [deg]");
                v.overall->GetYaxis()->SetTitle("#phi [deg]");
            }
            std::string out = outDir + "/" + v.saveName + "_" + selecteddetector + "_overall.png";
            c->SaveAs(out.c_str());
            delete c;
            std::cout << "Saved: " << out << std::endl;
        }

        for (size_t ti = 0; ti <= thetaCuts.size(); ++ti) {
            TCanvas *c = new TCanvas(("c_" + v.saveName + Form("_T%zu", ti)).c_str(), "", 1600, 1200);
            v.binHists[ti]->SetTitle((selecteddetector + " " + v.title + " in " + GetThetaBinLabel(ti, thetaCuts)).c_str());
            v.binHists[ti]->Draw("COLZ");
            if (v.name == "p:theta"){
                v.binHists[ti]->GetXaxis()->SetTitle("p [GeV]");
                v.binHists[ti]->GetYaxis()->SetTitle("#theta [deg]");
            } else if (v.name == "theta:p") {
                v.binHists[ti]->GetXaxis()->SetTitle("#theta [deg]");
                v.binHists[ti]->GetYaxis()->SetTitle("p [GeV]");
            } else if (v.name == "phi:theta") {
                v.binHists[ti]->GetXaxis()->SetTitle("#phi [deg]");
                v.binHists[ti]->GetYaxis()->SetTitle("#theta [deg]");
            } else if (v.name == "theta:phi") {
                v.binHists[ti]->GetXaxis()->SetTitle("#theta [deg]");
                v.binHists[ti]->GetYaxis()->SetTitle("#phi [deg]");
            }
            std::string out = outDir + "/" + v.saveName + "_" + selecteddetector + Form("_T%zu.png", ti);
            c->SaveAs(out.c_str());
            delete c;
            std::cout << "Saved: " << out << std::endl;
        }
    }

    timer.Stop();
    std::cout << "Time for Draw2DParticleKinematicsByThetaBins: ";
    timer.Print();
}

void DrawMCParticleKinematicsByThetaBins(
    const int &selectedPid,
    const std::vector<float> &thetaCuts,                                     // in degree
    const std::vector<std::tuple<std::string,std::string,std::string,int,double,double>> &plotVars,     // {name, nbins, xmin, xmax}
    const std::string &filename,
    const std::string &treename)
    {
    TStopwatch timer;
    timer.Start();

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df(treename, filename);
    //gStyle->SetOptStat(1110);

    struct VarInfo {
        std::string saveName;
        std::string title;  // for future use
        std::string name;
        int nbins;
        double xmin, xmax;
        TH1D* overall;
        std::vector<TH1D*> binHists;   // size = thetaCuts.size()+1
    };

    //==== create histograms ====
    std::vector<VarInfo> vars;
    for (auto &cfg : plotVars) {
        VarInfo v;
        v.saveName = std::get<0>(cfg);
        v.title    = std::get<1>(cfg);
        v.name  = std::get<2>(cfg);
        v.nbins = std::get<3>(cfg);
        v.xmin  = std::get<4>(cfg);
        v.xmax  = std::get<5>(cfg);
        v.overall = new TH1D((v.name+"_overall").c_str(), "", v.nbins, v.xmin, v.xmax);
        v.overall->SetDirectory(nullptr);
        v.binHists.resize(thetaCuts.size()+1);
        for (size_t ti=0; ti<=thetaCuts.size(); ++ti) {
            std::string hname = v.name + Form("_T%zu", ti);
            v.binHists[ti] = new TH1D(hname.c_str(), "", v.nbins, v.xmin, v.xmax);
            v.binHists[ti]->SetDirectory(nullptr);
        }
        vars.push_back(v);
    }

    //==== fill histograms ====
    std::vector<std::string> colNames = {"MC_Lund_pid","MC_Lund_px","MC_Lund_py","MC_Lund_pz"};

    df.Foreach([&](const RVec<int> &pid,
                   const RVec<float> &px,
                   const RVec<float> &py,
                   const RVec<float> &pz)  // passFid is not used but included for consistency
    {
        for (size_t i=0;i<pid.size();++i) {
            if (pid[i]!=selectedPid) continue;
            float p = std::sqrt(px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i]);
            float phi = std::atan2(py[i], px[i]);
            if (phi < 0) phi += 2 * M_PI;
            float theta = std::atan2(std::sqrt(px[i]*px[i] + py[i]*py[i]), pz[i]);
            float thetaDeg = theta*180.0/M_PI;
            int ti = GetThetaRegionIndex(thetaDeg, thetaCuts);
            for (auto &v: vars) {
                double value = NAN;
                if (v.name=="theta") value = theta*180.0/M_PI;
                else if (v.name=="phi") value = phi*180.0/M_PI; // ensure phi is in [0, 360)
                else if (v.name=="p") value = p;
                if (std::isnan(value)) continue;
                v.binHists[ti]->Fill(value);
                v.overall->Fill(value);
            }
        }
    }, colNames);

    //==== output directory ====
    std::string outDir = "MCParticleKinematicPlots";
    gSystem->Exec(("mkdir -p "+outDir).c_str());

    //==== draw and save ====
    for (auto &v: vars) {
        // overall
        {
            TCanvas *c = new TCanvas(("c_"+v.name+"_overall").c_str(),"",1600,1200);
            v.overall->SetTitle((v.title + " " + v.name).c_str());
            if (v.name == "phi" || v.name == "theta") v.overall->GetXaxis()->SetTitle((v.name+" [deg]").c_str());
            else v.overall->GetXaxis()->SetTitle((v.name+" GeV").c_str());
            v.overall->GetYaxis()->SetTitle("Counts");
            v.overall->Draw();
            std::string out = outDir + "/" + v.saveName + "_" + v.name + "_overall.png";
            c->SaveAs(out.c_str());
            delete c;
            std::cout << "Saved: " << out << std::endl;
        }
        // theta bins
        for (size_t ti=0; ti<=thetaCuts.size(); ++ti) {
            TCanvas *c = new TCanvas(("c_"+v.name+Form("_T%zu",ti)).c_str(),"",1600,1200);
            if (v.name == "phi" || v.name == "theta") v.binHists[ti]->GetXaxis()->SetTitle((v.name+ " [deg]").c_str());
            else v.binHists[ti]->GetXaxis()->SetTitle((v.name+" [GeV]").c_str()); 
            v.binHists[ti]->SetTitle((v.title+" "+v.name + " in " + GetThetaBinLabel(ti, thetaCuts)).c_str());
            v.binHists[ti]->GetYaxis()->SetTitle("Counts");
            v.binHists[ti]->Draw();
            std::string out = outDir + "/" + v.saveName + "_" + v.name + Form("_T%zu.png",ti);
            c->SaveAs(out.c_str());
            delete c;
            std::cout << "Saved: " << out << std::endl;
        }
    }

    timer.Stop();
    std::cout << "Time for DrawParticleKinematicsByThetaBins: ";
    timer.Print();
}

void DrawDeltaPByThetaBins(
    const int &selectedPid,
    const std::vector<float> &thetaCuts,
    const std::string selecteddetector,
    const std::vector<std::tuple<std::string,std::string,std::string,int,double,double,int,double,double>> &plotVars,
    const std::string &filename,
    const std::string &treename) {
    TStopwatch timer;
    timer.Start();

    ROOT::EnableImplicitMT();

    std::vector<double> thetaMidVec, aVec, bVec, cVec;
    std::vector<double> aErrVec, bErrVec, cErrVec;

    auto GetParticleName = [](int pid) -> std::string {
        if (pid == 11) return "electron";
        if (pid == 22) return "photon";
        if (pid == 2212) return "proton";
        return "pid" + std::to_string(pid);
    };
    std::string prefix = GetParticleName(selectedPid) + "_" + selecteddetector;


    ROOT::RDataFrame df(treename, filename);
    gStyle->SetOptStat(0);

    struct Var2DInfo {
        std::string saveName;
        std::string title;
        std::string name; // format: x:y (e.g., "deltaP:p")
        int nxbins;
        double xmin, xmax;
        int nybins;
        double ymin, ymax;
        TH2D* overall;
        std::vector<TH2D*> binHists;
    };

    std::vector<Var2DInfo> vars;
    for (auto &cfg : plotVars) {
        Var2DInfo v;
        v.saveName = std::get<0>(cfg);
        v.title    = std::get<1>(cfg);
        v.name     = std::get<2>(cfg);
        v.nxbins   = std::get<3>(cfg);
        v.xmin     = std::get<4>(cfg);
        v.xmax     = std::get<5>(cfg);
        v.nybins   = std::get<6>(cfg);
        v.ymin     = std::get<7>(cfg);
        v.ymax     = std::get<8>(cfg);

        v.overall = new TH2D((v.saveName+"_overall").c_str(), "", v.nxbins, v.xmin, v.xmax, v.nybins, v.ymin, v.ymax);
        v.overall->SetDirectory(nullptr);

        v.binHists.resize(thetaCuts.size() + 1);
        for (size_t ti = 0; ti <= thetaCuts.size(); ++ti) {
            std::string hname = v.saveName + Form("_T%zu", ti);
            v.binHists[ti] = new TH2D(hname.c_str(), "", v.nxbins, v.xmin, v.xmax, v.nybins, v.ymin, v.ymax);
            v.binHists[ti]->SetDirectory(nullptr);
        }
        vars.push_back(v);
    }

    std::vector<std::string> colNames = {"REC_Particle_pid", "REC_Particle_theta", "REC_Particle_phi", "REC_Particle_p", "REC_Particle_status", "REC_Particle_pass",
                                         "MC_Particle_pid", "MC_Particle_px", "MC_Particle_py", "MC_Particle_pz"};

    df.Foreach([&](const RVec<int> &pid,
                   const RVec<float> &theta,
                   const RVec<float> &phi,
                   const RVec<float> &p,
                   const RVec<short> &status,
                   const RVec<bool> &passPar,
                   const RVec<int> &mcpid,
                   const RVec<float> &mcpx,
                   const RVec<float> &mcpy,
                   const RVec<float> &mcpz)
    {
        std::vector<double> thetaCenters;
        std::vector<double> fitParamA, fitParamB, fitParamC;

        for (size_t i = 0; i < pid.size(); ++i) {
            if (pid[i] != selectedPid) continue;
            if (!passPar[i]) continue;
            if (GetDetectorPart(status[i]) != selecteddetector && selecteddetector != "ALL") continue;

            float thetaDeg = theta[i] * 180.0 / M_PI;
            int ti = GetThetaRegionIndex(thetaDeg, thetaCuts);

            float p_rec = p[i];
            float p_mc = NAN;
            for (size_t j = 0; j < mcpid.size(); ++j) {
                if (mcpid[j] == selectedPid) {
                    p_mc = std::sqrt(mcpx[j]*mcpx[j] + mcpy[j]*mcpy[j] + mcpz[j]*mcpz[j]);
                    break;
                }
            }
            if (std::isnan(p_mc)) continue;
            float deltaP = p_mc - p_rec;

            for (auto &v : vars) {
                double xval = NAN, yval = NAN;
                if (v.name == "deltaP:p") { xval = p_rec; yval = deltaP; }
                else if (v.name == "p:deltaP") { xval = deltaP; yval = p_rec; }
                else continue;

                if (std::isnan(xval) || std::isnan(yval)) continue;
                v.binHists[ti]->Fill(xval, yval);
                v.overall->Fill(xval, yval);
            }
        }
    }, colNames);

    std::string outDir = "ParticleDeltaPPlots";
    gSystem->Exec(("mkdir -p " + outDir).c_str());

    for (auto &v : vars) {
        {
            TCanvas *c = new TCanvas(("c_" + v.saveName + "_overall").c_str(), "", 3000, 1500);
            v.overall->SetTitle((selecteddetector + " " + v.title).c_str());
            gPad->SetLogz();
            v.overall->Draw("COLZ");
            //TGraph* gPeak = MakePeakGraph(v.overall,outDir + "/" + v.saveName + "_" + selecteddetector + "_overall.txt");
            //gPeak->Draw("PEZ SAME");
            //drawHist(v.overall, selecteddetector + " " + v.title, "p [GeV]", "#Delta p [GeV]");
            if (v.name == "deltaP:p") {
                v.overall->GetXaxis()->SetTitle("p [GeV]");
                v.overall->GetYaxis()->SetTitle("#Delta p [GeV]");
            } else if (v.name == "p:deltaP") {
                v.overall->GetXaxis()->SetTitle("#Delta p [GeV]");
                v.overall->GetYaxis()->SetTitle("p [GeV]");
            }
            std::string out = outDir + "/" + v.saveName + "_" + selecteddetector + "_overall.png";
            c->SaveAs(out.c_str());
            delete c;
            std::cout << "Saved: " << out << std::endl;
        }

        for (size_t ti = 0; ti <= thetaCuts.size(); ++ti) {
            TCanvas *c = new TCanvas(("c_" + v.saveName + Form("_T%zu", ti)).c_str(), "", 3000, 1500);
            v.binHists[ti]->SetTitle((selecteddetector + " " + v.title + " in " + GetThetaBinLabel(ti, thetaCuts)).c_str());
            gPad->SetLogz();
            v.binHists[ti]->Draw("COLZ");
            TGraph* gPeak2 = MakePeakGraph(v.binHists[ti],outDir + "/" + v.saveName + "_" + selecteddetector + Form("_T%zu.txt", ti));
            gPeak2->Draw("PEZ SAME");
  
            TF1* fitFunc = nullptr;
            std::string fitExpr;

            if (selectedPid == 2212) { // Proton
                if (selecteddetector == "FD") {
                    fitExpr = "[0] + [1]/x + [2]/(x*x)";
                } else if (selecteddetector == "CD") {
                    fitExpr = "[0] + [1]*x + [2]*x*x";  // a + bx + cx^2
                } else {
                    std::cerr << "Unknown detector for proton: " << selecteddetector << "\n";
                    continue;
                }
            } else if (selectedPid == 11 || selectedPid == 22) { // Electron or Photon
                fitExpr = "[0] + [1]*x + [2]*x*x"; 
            } else {
                std::cerr << "Unknown pid for fitting: " << selectedPid << "\n";
                fitExpr = "[0] + [1]*x + [2]*x*x"; // fallback
            }

            fitFunc = new TF1("fitFunc", fitExpr.c_str(), v.binHists[ti]->GetXaxis()->GetXmin(), v.binHists[ti]->GetXaxis()->GetXmax());


            if (gPeak2->GetN() > 2) {
                gPeak2->Fit(fitFunc, "Q");  // Quiet mode

                std::ofstream foutFit(outDir + "/" + v.saveName + "_" + selecteddetector + Form("_T%zu_fit.txt", ti));
                foutFit << "# Fit expression: " << fitExpr << "\n";
                foutFit << "# Fit parameters:\n";
                for (int i = 0; i < fitFunc->GetNpar(); ++i)
                    foutFit << "p" << i << " = " << fitFunc->GetParameter(i)
                            << " ± " << fitFunc->GetParError(i) << "\n";
                foutFit.close();

                fitFunc->SetLineColor(kBlue + 2);
                fitFunc->SetLineWidth(2);
                fitFunc->Draw("SAME");

                // === 记录每个 bin 的 theta center 和 拟合参数 ===
                double thetaMid;
                if (ti == 0)
                    thetaMid = thetaCuts[0];
                else if (ti == thetaCuts.size())
                    thetaMid = thetaCuts.back();
                else
                    thetaMid = 0.5 * (thetaCuts[ti] + thetaCuts[ti - 1]);

                if (fitFunc->GetNpar() >= 3 && ti < thetaCuts.size()-1) {
                    thetaMidVec.push_back(thetaMid);
                    aVec.push_back(fitFunc->GetParameter(0));
                    bVec.push_back(fitFunc->GetParameter(1));
                    cVec.push_back(fitFunc->GetParameter(2));
                    aErrVec.push_back(fitFunc->GetParError(0));
                    bErrVec.push_back(fitFunc->GetParError(1));
                    cErrVec.push_back(fitFunc->GetParError(2));

                }
            }

            //drawHist(v.binHists[ti],selecteddetector + " " + v.title + " in " + GetThetaBinLabel(ti, thetaCuts),"p [GeV]", "#Delta p [GeV]");

            if (v.name == "deltaP:p") {
                v.binHists[ti]->GetXaxis()->SetTitle("p [GeV]");
                v.binHists[ti]->GetYaxis()->SetTitle("#Delta p [GeV]");
            } else if (v.name == "p:deltaP") {
                v.binHists[ti]->GetXaxis()->SetTitle("#Delta p [GeV]");
                v.binHists[ti]->GetYaxis()->SetTitle("p [GeV]");
            }
            std::string out = outDir + "/" + v.saveName + "_" + selecteddetector + Form("_T%zu.png", ti);
            c->SaveAs(out.c_str());
            delete c;
            std::cout << "Saved: " << out << std::endl;
        }
    }
    std::string fitParamOutDir = outDir + "/ParamFits";
    gSystem->Exec(("mkdir -p " + fitParamOutDir).c_str());

    auto PlotParamVsTheta = [&](const std::vector<double>& theta,
                             const std::vector<double>& param,
                             const std::vector<double>& errors,
                             const std::string& pname,
                             const std::string& fitExpr) {
        if (theta.size() < 3) return;

        TCanvas *c = new TCanvas(("c_" + prefix + "_" + pname + "_vs_theta").c_str(), "", 1800, 1200);
        TGraphErrors* g = new TGraphErrors(theta.size());
        for (size_t i = 0; i < theta.size(); ++i) {
            g->SetPoint(i, theta[i], param[i]);
            g->SetPointError(i, 0.0, errors[i]); // X误差设为0，Y误差来自拟合
        }

        TF1* f = new TF1(("fit_" + pname).c_str(), fitExpr.c_str(),
                     *std::min_element(theta.begin(), theta.end()),
                     *std::max_element(theta.begin(), theta.end()));
        g->Fit(f, "Q");

        g->SetTitle((prefix + ": " + pname + " vs #theta").c_str());
        g->GetXaxis()->SetTitle("#theta [deg]");
        g->GetYaxis()->SetTitle((pname + " value").c_str());
        g->SetMarkerStyle(20);
        g->SetMarkerSize(1.2);
        g->Draw("AP");

        f->SetLineColor(kRed + 1);
        f->SetLineWidth(2);
        f->Draw("SAME");
        g->GetYaxis()->SetRangeUser(-0.05, 0.05);


        std::string imgPath = fitParamOutDir + "/" + prefix + "_" + pname + "_vs_theta.png";
        c->SaveAs(imgPath.c_str());
        delete c;
        std::cout << "Saved: " << imgPath << std::endl;

        std::ofstream fout(fitParamOutDir + "/" + prefix + "_" + pname + "_vs_theta.txt");
        fout << "# Fit expression: " << fitExpr << "\n";
        for (int i = 0; i < f->GetNpar(); ++i)
        fout << "p" << i << " = " << f->GetParameter(i)
             << " ± " << f->GetParError(i) << "\n";
        fout.close();
    };

    if (selecteddetector == "FD") {
        PlotParamVsTheta(thetaMidVec, aVec, aErrVec, "A_p", "[0] + [1]*x");
        PlotParamVsTheta(thetaMidVec, bVec, bErrVec, "B_p", "[0] + [1]*x");
        PlotParamVsTheta(thetaMidVec, cVec, cErrVec, "C_p", "[0] + [1]*x");
    } else if (selecteddetector == "CD") {
        PlotParamVsTheta(thetaMidVec, aVec, aErrVec, "A_p", "[0] + [1]*x");
        PlotParamVsTheta(thetaMidVec, bVec, bErrVec, "B_p", "[0] + [1]*x");
        PlotParamVsTheta(thetaMidVec, cVec, cErrVec, "C_p", "[0] + [1]*x");
    }



    timer.Stop();
    std::cout << "Time for DrawDeltaPByThetaBins: ";
    timer.Print();
}




//================ example driver =================
void analysisMomentumCorrection() {
    //std::string path = "../build/rgk7546dvcsmcAll/";
    std::string path = "/w/hallb-scshelf2102/clas12/yijie/clas12ana/analysis316/DISANA/build/rgk7546clasdismc/";
    std::string filename = path + "dfSelected_afterFid.root";
    std::string filenameCorrected = path + "dfSelected_afterFid_afterCorr.root";
    std::string treename = "dfSelected_afterFid";
    std::string treenameCorrected = "dfSelected_afterFid_afterCorr";

    std::vector<float> thetaCutsFDelectron = {10,15,20,25};
    std::vector<float> thetaCutsFDproton = {13,14,15,16,17,18,19,20,21,22,23,24,25,
                                            26,27,28,29,30,31,32,33,34,35,36,37,38,39,40};
    std::vector<float> thetaCutsFDphoton = {10,15,20,25};

    std::vector<float> thetaCutsCDproton = {40,42,44,46,48,50,52,54,56,58,60};

    std::vector<float> thetaCutsFTelectron = {3,3.5,4};
    std::vector<float> thetaCutsFTphoton = {3,3.5,4};

    std::vector<std::tuple<std::string,std::string,std::string,int,double,double>> plotVarsFDelectron = {
        {"electron","electron","theta", 500, 0, 60},
        {"electron","electron","phi",   500, 0, 360},
        {"electron","electron","p",     500, 0, 8.0}
    };
    std::vector<std::tuple<std::string,std::string,std::string,int,double,double>> plotVarsFDproton = {
        {"proton","proton","theta", 500, 0, 60},
        {"proton","proton","phi",   500, 0, 360},
        {"proton","proton","p",     500, 0, 3.0}
    };
    std::vector<std::tuple<std::string,std::string,std::string,int,double,double>> plotVarsFDphoton = {
        {"photon","photon","theta", 500, 0, 60},
        {"photon","photon","phi",   500, 0, 360},
        {"photon","photon","p",     500, 0, 8.0}
    };
    std::vector<std::tuple<std::string,std::string,std::string,int,double,double>> plotVarsCDproton = {
        {"proton","proton","theta", 500, 20, 150},
        {"proton","proton","phi",   500, 0, 360},
        {"proton","proton","p",     500, 0, 2.0}
    };
    std::vector<std::tuple<std::string,std::string,std::string,int,double,double>> plotVarsFTelectron = {
        {"electron","electron","theta", 500, 0, 10},
        {"electron","electron","phi",   500, 0, 360},
        {"electron","electron","p",     500, 0, 8.0}
    };
    std::vector<std::tuple<std::string,std::string,std::string,int,double,double>> plotVarsFTphoton = {
        {"photon","photon","theta", 500, 0, 10},
        {"photon","photon","phi",   500, 0, 360},
        {"photon","photon","p",     500, 0, 8.0}
    };
    std::vector<std::tuple<std::string,std::string,std::string,int,double,double>> plotVarsAllelectron = {
        {"electron","electron","theta", 500, 0, 60},
        {"electron","electron","phi",   500, 0, 360},
        {"electron","electron","p",     500, 0, 8.0},
        {"electron","electron","beta",  500, 0, 2.0}
    };
    std::vector<std::tuple<std::string,std::string,std::string,int,double,double>> plotVarsAllphoton = {
        {"photon","photon","theta", 500, 0, 60},
        {"photon","photon","phi",   500, 0, 360},
        {"photon","photon","p",     500, 0, 8.0},
        {"photon","photon","beta",  500, 0, 2.0}
    };
    std::vector<std::tuple<std::string,std::string,std::string,int,double,double>> plotVarsAllproton = {
        {"proton","proton","theta", 500, 0, 150},
        {"proton","proton","phi",   500, 0, 360},
        {"proton","proton","p",     500, 0, 3.0},
        {"proton","proton","beta",  500, 0, 2.0}
    };
/*
    DrawParticleKinematicsByThetaBins(22, {0}, "ALL", plotVarsAllphoton, filename, treename); // photon
    DrawParticleKinematicsByThetaBins(11, {0}, "ALL", plotVarsAllelectron, filenameCorrected, treenameCorrected); // electron
    DrawParticleKinematicsByThetaBins(2212, {0}, "ALL", plotVarsAllproton, filenameCorrected, treenameCorrected); // proton
*/
/*

    DrawParticleKinematicsByThetaBins(11, thetaCutsFDelectron, "FD", plotVarsFDelectron, filename, treename); // electron
    DrawParticleKinematicsByThetaBins(2212, thetaCutsFDproton, "FD", plotVarsFDproton, filename, treename); // proton
    DrawParticleKinematicsByThetaBins(22, thetaCutsFDphoton, "FD", plotVarsFDphoton, filename, treename); // photon

    DrawParticleKinematicsByThetaBins(2212, thetaCutsCDproton, "CD", plotVarsCDproton, filename, treename); // proton

    DrawParticleKinematicsByThetaBins(11, thetaCutsFTelectron, "FT", plotVarsFTelectron, filename, treename); // electron
    DrawParticleKinematicsByThetaBins(22, thetaCutsFTphoton, "FT", plotVarsFTphoton, filename, treename); // photon

    DrawParticleKinematicsByThetaBins(11, {0}, "ALL", plotVarsAllelectron, filename, treename); // electron
    DrawParticleKinematicsByThetaBins(22, {0}, "ALL", plotVarsAllphoton, filename, treename); // photon
    DrawParticleKinematicsByThetaBins(2212, {0}, "ALL", plotVarsAllproton, filename, treename); // proton

    DrawMCParticleKinematicsByThetaBins(11, {0}, plotVarsAllelectron, filename, treename); // electron
    DrawMCParticleKinematicsByThetaBins(2212, {0}, plotVarsAllproton, filename, treename); // proton
    DrawMCParticleKinematicsByThetaBins(22, {0}, plotVarsAllphoton, filename, treename); // photon

    Draw2DParticleKinematicsByThetaBins(11,{0},"FD",{{"electron_p_vs_theta","electron p vs theta","p:theta",500,0,10,500,0,40}},filename,treename);
    Draw2DParticleKinematicsByThetaBins(11,{0},"FD",{{"electron_phi_vs_theta","electron phi vs theta","phi:theta",500,0,360,500,0,40}},filename,treename);
    Draw2DParticleKinematicsByThetaBins(11,{0},"FT",{{"electron_p_vs_theta","electron p vs theta","p:theta",500,0,10,500,0,10}},filename,treename);
    Draw2DParticleKinematicsByThetaBins(11,{0},"FT",{{"electron_phi_vs_theta","electron phi vs theta","phi:theta",500,0,360,500,0,10}},filename,treename);
  Draw2DParticleKinematicsByThetaBins(11,{0},"ALL",{{"electron_p_vs_theta","electron p vs theta","p:theta",500,0,10,500,0,40}},filename,treename);
*/
/*
    Draw2DParticleKinematicsByThetaBins(2212,{0},"FD",{{"proton_p_vs_theta","proton p vs theta","p:theta",500,0,3,500,0,60}},filename,treename);
    Draw2DParticleKinematicsByThetaBins(2212,{0},"FD",{{"proton_phi_vs_theta","proton phi vs theta","phi:theta",500,0,360,500,0,60}},filename,treename);
    Draw2DParticleKinematicsByThetaBins(2212,{0},"CD",{{"proton_p_vs_theta","proton p vs theta","p:theta",500,0,2,500,20,150}},filename,treename);
    Draw2DParticleKinematicsByThetaBins(2212,{0},"CD",{{"proton_phi_vs_theta","proton phi vs theta","phi:theta",500,0,360,500,20,150}},filename,treename);
    Draw2DParticleKinematicsByThetaBins(2212,{0},"ALL",{{"proton_p_vs_theta","proton p vs theta","p:theta",500,0,3,500,0,150}},filename,treename);
  */
/*
    Draw2DParticleKinematicsByThetaBins(22,{0},"FD",{{"photon_p_vs_theta","photon p vs theta","p:theta",500,0,8,500,0,40}},filenameCorrected,treenameCorrected);
    Draw2DParticleKinematicsByThetaBins(22,{0},"FD",{{"photon_phi_vs_theta","photon phi vs theta","phi:theta",500,0,360,500,0,40}},filenameCorrected,treenameCorrected);
    Draw2DParticleKinematicsByThetaBins(22,{0},"FT",{{"photon_p_vs_theta","photon p vs theta","p:theta",500,0,8,500,0,10}},filenameCorrected,treenameCorrected);
    Draw2DParticleKinematicsByThetaBins(22,{0},"FT",{{"photon_phi_vs_theta","photon phi vs theta","phi:theta",500,0,360,500,0,10}},filenameCorrected,treenameCorrected);
    Draw2DParticleKinematicsByThetaBins(22,{0},"ALL",{{"photon_p_vs_theta","photon p vs theta","p:theta",500,0,8,500,0,40}},filenameCorrected,treenameCorrected);
    Draw2DParticleKinematicsByThetaBins(22,{0},"ALL",{{"photon_phi_vs_theta","photon phi vs theta","phi:theta",500,0,360,500,0,40}},filenameCorrected,treenameCorrected);
*/
/*
    DrawDeltaPByThetaBins(11,thetaCutsFDelectron,"FD",{{"electron_deltaP_vs_p","electron #Delta p vs p","deltaP:p",500,0,8,500,-0.03,0.03}},filename,treename);
    DrawDeltaPByThetaBins(11,thetaCutsFTelectron,"FT",{{"electron_deltaP_vs_p","electron #Delta p vs p","deltaP:p",500,0,8,500,-1.0,1.0}},filename,treename);
    DrawDeltaPByThetaBins(11,{0},"ALL",{{"electron_deltaP_vs_p","electron #Delta p vs p","deltaP:p",500,0,8,500,-0.1,0.1}},filename,treename);
  */
   
    DrawDeltaPByThetaBins(2212,thetaCutsFDproton,"FD",{{"proton_deltaP_vs_p","proton #Delta p vs p","deltaP:p",100,0,2,100,-0.1,0.1}},filename,treename);
    DrawDeltaPByThetaBins(2212,thetaCutsCDproton,"CD",{{"proton_deltaP_vs_p","proton #Delta p vs p","deltaP:p",100,0,2,100,-0.2,0.2}},filename,treename);
    DrawDeltaPByThetaBins(2212,{0},"ALL",{{"proton_deltaP_vs_p","proton #Delta p vs p","deltaP:p",100,0,2,100,-0.1,0.1}},filename,treename);
    //DrawDeltaPByThetaBins(2212,thetaCutsFDproton,"FD",{{"proton_deltaP_vs_pcorr","corrected proton #Delta p vs p","deltaP:p",100,0,2,100,-0.1,0.1}},filenameCorrected,treenameCorrected);
    //DrawDeltaPByThetaBins(2212,thetaCutsCDproton,"CD",{{"proton_deltaP_vs_pcorr","corrected proton #Delta p vs p","deltaP:p",100,0,2,100,-0.2,0.2}},filenameCorrected,treenameCorrected);

    /*
    DrawDeltaPByThetaBins(22,thetaCutsFDphoton,"FD",{{"photon_deltaP_vs_p","photon #Delta p vs p","deltaP:p",500,0,8,500,-0.5,0.5}},filename,treename);
    DrawDeltaPByThetaBins(22,thetaCutsFTphoton,"FT",{{"photon_deltaP_vs_p","photon #Delta p vs p","deltaP:p",500,0,8,500,-0.2,0.2}},filename,treename);
    DrawDeltaPByThetaBins(22,{0},"ALL",{{"photon_deltaP_vs_p","photon #Delta p vs p","deltaP:p",500,0,8,500,-0.5,0.5}},filename,treename);
*/


    gApplication->Terminate(0);
}
