#ifndef COLUMNS_H_
#define COLUMNS_H_

#include <vector>
#include <string>
#include <algorithm>
#include <limits>
#include <type_traits>

class Columns {
 public:
  // For two-column logical AND
  static std::function<std::vector<int>(const std::vector<int>&, const std::vector<int>&)> LogicalAND2();
};

// CombineColumns 函数声明
template <typename... Args>
std::vector<std::string> CombineColumns(const Args&... vectors);
#endif // COLUMNS_H_