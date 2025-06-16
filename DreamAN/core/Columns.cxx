#include "Columns.h"

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
