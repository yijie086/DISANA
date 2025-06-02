#include "Columns.h"

// CombineColumns 函数定义
template <typename... Args>
std::vector<std::string> CombineColumns(const Args&... vectors) {
    std::vector<std::string> combined;
    // 使用 fold expression 合并所有输入的 std::vector<std::string>
    (combined.insert(combined.end(), vectors.begin(), vectors.end()), ...);
    return combined;
}

// 显式实例化模板
template std::vector<std::string> CombineColumns(const std::vector<std::string>&);
template std::vector<std::string> CombineColumns(const std::vector<std::string>&, const std::vector<std::string>&);
template std::vector<std::string> CombineColumns(const std::vector<std::string>&, const std::vector<std::string>&, const std::vector<std::string>&);