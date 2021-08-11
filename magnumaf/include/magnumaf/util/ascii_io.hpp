#pragma once
#include "arrayfire.h"
#include "mesh.hpp"
#include <string>
#include <utility>

namespace magnumaf {

void write_ascii(const af::array& a, const Mesh& mesh, const std::string& filename, bool verbose = true,
                 int precision = 18);
std::pair<af::array, Mesh> read_ascii(const std::string& filename, bool verbose = true);

} // namespace magnumaf
