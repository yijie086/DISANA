#ifndef COLUMNS_H_
#define COLUMNS_H_

#include <vector>
#include <string>
#include <algorithm>
#include <limits>
#include <type_traits>
#include <functional>

class Columns {
 public:
  // For two-column logical AND
  static std::function<std::vector<int>(const std::vector<int>&, const std::vector<int>&)> LogicalAND2();
};

// CombineColumns: Merges multiple std::vector<std::string> into one
template <typename... Args>
std::vector<std::string> CombineColumns(const Args&... vectors) {
  std::vector<std::string> combined;
  (combined.insert(combined.end(), vectors.begin(), vectors.end()), ...);
  return combined;
}

#endif // COLUMNS_H_
