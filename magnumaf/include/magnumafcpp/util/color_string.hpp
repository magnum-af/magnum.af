#pragma once
#include <string>

//! Colorizing strings primarily for std::cout.
namespace magnumaf::color_string {

///> Red is for warnings.
inline std::string red(const std::string& str) { return "\033[;31m" + str + "\033[0m"; }

///> Bold red is for errors.
inline std::string bold_red(const std::string& str) { return "\033[1;31m" + str + "\033[0m"; }

///> Green is for regular infos.
inline std::string bold_green(const std::string& str) { return "\033[1;32m" + str + "\033[0m"; }

///> Bold green is for sucess infos.
inline std::string green(const std::string& str) { return "\033[;32m" + str + "\033[0m"; }

///> Printing a colorized Info string
inline const char* info() { return "\33[0;32mInfo:\33[0m"; }

///> Printing a colorized Warning string
inline const char* warning() { return "\33[1;31mWarning:\33[0m"; }

} // namespace magnumaf::color_string
