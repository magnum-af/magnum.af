#pragma once
#include "arrayfire.h"
#include <iostream>
#include <ostream>

namespace magnumafcpp {

struct NonequispacedMesh {
    NonequispacedMesh(unsigned nx, unsigned ny, double dx, double dy, std::vector<double> z_spacing);
    unsigned nx, ny, nz;               //!< Number of cells in x, y, z
    double dx, dy;                     //!< Distance between equidistant x, y cells
    std::vector<double> z_spacing;     //
};

std::ostream& operator<<(std::ostream& os, const NonequispacedMesh& nemesh);

inline af::dim4 dims_vector(NonequispacedMesh nemesh) { return af::dim4(nemesh.nx, nemesh.ny, nemesh.nz, 3); }

} // namespace magnumafcpp
