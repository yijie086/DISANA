#include <TFile.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <iostream>
#include <string>
#include <vector>

void DrawProfiles() {
    // 打开 ROOT 文件
    TFile* file = TFile::Open("../build/test.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open the ROOT file!" << std::endl;
        return;
    }
    std::cout << "Successfully opened the ROOT file." << std::endl;

    // 遍历 R1 到 R3
    for (int region = 1; region <= 3; ++region) {
        std::string region_str = "R" + std::to_string(region);

        // 遍历 sector 1 到 sector 6
        for (int sector = 1; sector <= 6; ++sector) {
            // 定义图的名字
            std::vector<std::string> profile_names = {
                "theta5to10_sector" + std::to_string(sector) + "_" + region_str + "_chi2perndf_vs_edge_profile",
                "theta10to15_sector" + std::to_string(sector) + "_" + region_str + "_chi2perndf_vs_edge_profile",
                "theta15to20_sector" + std::to_string(sector) + "_" + region_str + "_chi2perndf_vs_edge_profile",
                "theta20to25_sector" + std::to_string(sector) + "_" + region_str + "_chi2perndf_vs_edge_profile",
                "theta25to40_sector" + std::to_string(sector) + "_" + region_str + "_chi2perndf_vs_edge_profile"
            };

            // 定义图例标签
            std::vector<std::string> legend_labels = {
                "5-10 degrees",
                "10-15 degrees",
                "15-20 degrees",
                "20-25 degrees",
                "25-40 degrees"
            };

            // 创建画布
            TCanvas* canvas = new TCanvas(("canvas_" + region_str + "_sector" + std::to_string(sector)).c_str(),
                                          ("Chi2/NDF vs Edge " + region_str + " Sector " + std::to_string(sector)).c_str(),
                                          1920, 1080);
            canvas->SetGrid();

            // 创建图例
            TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9); // 图例位置
            legend->SetBorderSize(0); // 无边框
            legend->SetFillStyle(0);  // 透明背景

            // 全局关闭统计框
            gStyle->SetOptStat(0);

            // 颜色数组
            int colors[] = {kRed, kBlue, kGreen, kMagenta, kCyan, kOrange};

            // 读取并绘制每个 TProfile
            for (int i = 0; i < 5; ++i) {
                TProfile* profile = dynamic_cast<TProfile*>(file->Get(profile_names[i].c_str()));
                if (!profile) {
                    std::cerr << "Error: Cannot find profile " << profile_names[i] << " in the file!" << std::endl;
                    continue;
                }
                for (int bin = 1; bin <= profile->GetNbinsX(); ++bin) {
                    double x = profile->GetBinCenter(bin); // 获取 bin 的中心值
                    if (x > 15 && i == 0 && region == 1) {
                        profile->SetBinContent(bin, 0); // 将 bin 内容设置为 0
                        profile->SetBinError(bin, 0);  // 将 bin 的误差设置为 0
                    }
                }

                // 设置线条颜色和样式
                profile->SetLineColor(colors[i]);
                profile->SetLineWidth(2);
                profile->GetYaxis()->SetRangeUser(0, 200);
                profile->SetStats(0);
                profile->SetTitle(""); // 清除默认标题
                profile->GetXaxis()->SetTitle("edge/cm");       // 设置 x 轴标签
                profile->GetYaxis()->SetTitle("<Chi2/NDF>");    // 设置 y 轴标签

                // 绘制图像
                if (i == 0) {
                    profile->Draw(""); // 第一个图用 "HIST" 绘制
                } else {
                    profile->Draw("SAME"); // 其余图用 "HIST SAME" 绘制
                }

                // 添加到图例
                legend->AddEntry(profile, legend_labels[i].c_str(), "l");
            }

            // 绘制图例
            legend->Draw();

            // 添加标题
            TLatex latex;
            latex.SetTextSize(0.04); // 设置字体大小
            latex.SetTextAlign(22);  // 设置对齐方式（22 表示居中）
            latex.DrawLatexNDC(0.5, 0.95, ("electron Chi2/NDF vs Edge " + region_str + " sector " + std::to_string(sector)).c_str());

            // 保存画布为图片
            canvas->SaveAs(("chi2perndf_vs_edge_profiles_" + region_str + "_sector" + std::to_string(sector) + ".png").c_str());

            // 删除画布
            delete canvas;
        }
    }

    // 关闭文件
    file->Close();
    delete file;
}