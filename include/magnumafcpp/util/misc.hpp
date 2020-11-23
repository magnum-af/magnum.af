#pragma once
#include <string>
#include <sys/stat.h>

namespace magnumafcpp {

///> Colorizing a std::string for std::cout. Red is for warnings.
inline std::string red(const std::string str) { return "\033[;31m" + str + "\033[0m"; }
///> Colorizing a std::string for std::cout. Bold red is for errors.
inline std::string bold_red(const std::string str) { return "\033[1;31m" + str + "\033[0m"; }
///> Colorizing a std::string for std::cout. Green is for regular infos.
inline std::string bold_green(const std::string str) { return "\033[1;32m" + str + "\033[0m"; }
///> Colorizing a std::string for std::cout. Bold green is for sucess infos.
inline std::string green(const std::string str) { return "\033[;32m" + str + "\033[0m"; }

inline const char* Info(void) { return "\33[0;32mInfo:\33[0m"; }
inline const char* Warning(void) { return "\33[1;31mWarning:\33[0m"; }
} // namespace magnumafcpp
