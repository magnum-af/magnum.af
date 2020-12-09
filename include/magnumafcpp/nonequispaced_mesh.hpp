#pragma once
#include "arrayfire.h"
#include <cstddef>
#include <iostream>
#include <ostream>

namespace magnumafcpp {

struct NonequiMesh {
    NonequiMesh(std::size_t nx, std::size_t ny, double dx, double dy, std::vector<double> z_spacing)
        : nx(nx), ny(ny), nz(z_spacing.size()), dx(dx), dy(dy), z_spacing(z_spacing) {}
    std::size_t nx, ny, nz;        //!< Number of cells in x, y, z
    double dx, dy;                 //!< Distance between equidistant x, y cells
    std::vector<double> z_spacing; //
};

inline std::ostream& operator<<(std::ostream& os, const NonequiMesh& nemesh) {
    os << "nx=" << nemesh.nx << " ny=" << nemesh.ny << " nz=" << nemesh.nz << " dx=" << nemesh.dx << " dy=" << nemesh.dy
       << " dz: ";
    for (auto const& dz : nemesh.z_spacing) {
        os << dz << " ";
    }
    return os;
}

inline af::dim4 dims_vector(NonequiMesh nemesh) { return af::dim4(nemesh.nx, nemesh.ny, nemesh.nz, 3); }

} // namespace magnumafcpp
