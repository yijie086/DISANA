#include "Columns.h"

// CombineColumns 函数定义
template <typename... Args>
std::vector<std::string> CombineColumns(const Args&... vectors) {
  std::vector<std::string> combined;
  // 使用 fold expression 合并所有输入的 std::vector<std::string>
  (combined.insert(combined.end(), vectors.begin(), vectors.end()), ...);
  return combined;
}

std::function<std::vector<int>(const std::vector<int>&, const std::vector<int>&)>
Columns::LogicalAND2() {
  return [](const std::vector<int>& a, const std::vector<int>& b) {
    size_t N = std::min(a.size(), b.size());
    std::vector<int> result(N, 1);
    for (size_t i = 0; i < N; ++i)
      result[i] = a[i] && b[i];
    return result;
  };
}


// 显式实例化模板
template std::vector<std::string> CombineColumns(const std::vector<std::string>&);
template std::vector<std::string> CombineColumns(const std::vector<std::string>&, const std::vector<std::string>&);
template std::vector<std::string> CombineColumns(const std::vector<std::string>&, const std::vector<std::string>&, const std::vector<std::string>&);

template std::vector<std::string> CombineColumns(const std::vector<std::string>&, const std::vector<std::string>&, const std::vector<std::string>&, const std::vector<std::string>&, const std::vector<std::string>&, const std::vector<std::string>&);