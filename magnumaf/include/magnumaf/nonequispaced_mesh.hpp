#pragma once
#include "arrayfire.h"
#include <cstddef>
#include <iostream>
#include <ostream>

namespace magnumaf {

using std::size_t;

struct NonequiMesh {
    NonequiMesh(size_t nx, size_t ny, double dx, double dy, std::vector<double> z_spacing);
    size_t nx, ny;                 //!< Number of cells in x, y
    size_t nz;                     //!< Number of cells in z, determined by length of z_spacing
    double dx, dy;                 //!< Distance between equidistant x, y cells in [m]
    std::vector<double> z_spacing; //!< Thickness of each layer along z-axis in [m]
};

inline NonequiMesh::NonequiMesh(size_t nx_, size_t ny_, double dx_, double dy_, std::vector<double> z_spacing_)
    : nx(nx_), ny(ny_), nz(z_spacing_.size()), dx(dx_), dy(dy_), z_spacing(z_spacing_) {}

inline std::ostream& operator<<(std::ostream& os, const NonequiMesh& nemesh) {
    os << "nx=" << nemesh.nx << " ny=" << nemesh.ny << " nz=" << nemesh.nz << " dx=" << nemesh.dx << " dy=" << nemesh.dy
       << " dz: ";
    for (const auto& dz : nemesh.z_spacing) {
        os << dz << " ";
    }
    return os;
}

namespace nemesh {
/// Dimension for vector field on nemesh, i.e. af::dim4(nx, ny, nz, 3)
inline af::dim4 dims_s(NonequiMesh nemesh) { return af::dim4(nemesh.nx, nemesh.ny, nemesh.nz, 1); }
inline af::dim4 dims_v(NonequiMesh nemesh) { return af::dim4(nemesh.nx, nemesh.ny, nemesh.nz, 3); }
} // namespace nemesh

} // namespace magnumaf
