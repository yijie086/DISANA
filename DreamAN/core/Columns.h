#ifndef COLUMNS_H_
#define COLUMNS_H_

#include <vector>
#include <string>

// CombineColumns 函数声明
template <typename... Args>
std::vector<std::string> CombineColumns(const Args&... vectors);

#endif // COLUMNS_H_