#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <iostream>

void SplitRDF(const std::string& inputFile = "../build/qadbcut/dfSelected_afterFid_afterCorr.root", const std::string& treeName = "dfSelected_afterFid_afterCorr") {
  // 开启多线程（可选）
  ROOT::EnableImplicitMT();

  // 直接用完整命名空间
  ROOT::RDataFrame df(treeName, inputFile);

  // 按 run 分割
  auto df_low  = df.Filter("RUN_config_run <= 5757", "run <= 5757");
  auto df_high = df.Filter("RUN_config_run > 5757",  "run > 5757");

  // 打印统计信息
  std::cout << "Total entries: " << *df.Count() << std::endl;
  std::cout << "Entries run<=5757: " << *df_low.Count() << std::endl;
  std::cout << "Entries run>5757 : " << *df_high.Count() << std::endl;

  // 保存两个文件
  df_low.Snapshot(treeName, "output_runLE5757.root");
  df_high.Snapshot(treeName, "output_runGT5757.root");

  std::cout << "✅ Saved: output_runLE5757.root and output_runGT5757.root" << std::endl;
}
