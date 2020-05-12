#pragma once
#include "arrayfire.h"
#include "mesh.hpp"
#include <string>

namespace magnumafcpp {

void write_ascii(const af::array& a, const Mesh& mesh, std::string filename,
                 bool verbose = true, int precision = 18);
void read_ascii(af::array& a, Mesh& mesh, std::string filename,
                bool verbose = true);

} // namespace magnumafcpp