#pragma once
#include "arrayfire.h"
#include "mesh.hpp"
#include <string>
#include <utility>

namespace magnumafcpp {

void write_ascii(const af::array& a, const Mesh& mesh, std::string filename, bool verbose = true, int precision = 18);
std::pair<af::array, Mesh> read_ascii(std::string filename, bool verbose = true);

} // namespace magnumafcpp
