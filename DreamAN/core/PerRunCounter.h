// ===== PerRunCounter.h =====
#pragma once
#include <ROOT/RDataFrame.hxx>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <algorithm>
#include <type_traits>
#include <string>

// 统计 df 上指定 run 列的频次，并写出 CSV: run,count
// 返回 <run, count> 的升序向量。
// TRun 必须是整型（int, long, int64_t ...）
template <typename TRun>
std::vector<std::pair<TRun, long long>>
CountPerRunAndWriteCSV(ROOT::RDF::RNode df,
                       const std::string& run_col,
                       const std::string& out_csv_path)
{
  static_assert(std::is_integral<TRun>::value, "TRun must be an integral type");

  // 触发执行并取出 run 列
  auto runs_res = df.Take<TRun>(run_col);
  const std::vector<TRun>& runs = runs_res.GetValue();

  // 计数
  std::unordered_map<TRun, long long> cnt;
  cnt.reserve(runs.size()/8 + 8);
  for (auto r : runs) ++cnt[r];

  // 升序整理
  std::vector<std::pair<TRun, long long>> items(cnt.begin(), cnt.end());
  std::sort(items.begin(), items.end(),
            [](const auto& a, const auto& b){ return a.first < b.first; });

  // 写 CSV
  std::ofstream ofs(out_csv_path);
  ofs << "run,events_afterFid\n";
  for (const auto& kv : items) ofs << kv.first << "," << kv.second << "\n";
  ofs.close();

  return items; // 如不需要返回，可改为 void
}
