#pragma once
#include "arrayfire.h"
#include <iostream>

namespace magnumafcpp {

struct Mesh {
    Mesh(unsigned nx, unsigned ny, unsigned nz, double dx, double dy, double dz)
        : nx(nx), ny(ny), nz(nz), dx(dx), dy(dy), dz(dz) {}

    unsigned nx, ny, nz; ///< Number of cells in x, y, z
    double dx, dy, dz;   ///< Distance between cells
};

inline std::ostream& operator<<(std::ostream& os, const Mesh& mesh) {
    os << "nx=" << mesh.nx << " ny=" << mesh.ny << " nz=" << mesh.nz << " dx=" << mesh.dx << " dy=" << mesh.dy
       << " dz=" << mesh.dz;
    return os;
}

// Dimension for scalar field on mesh, i.e. [nx, ny, nz, 1]
inline af::dim4 dims_scalar(Mesh mesh) { return af::dim4(mesh.nx, mesh.ny, mesh.nz, 1); }

// Dimension for vector field on mesh, i.e. [nx, ny, nz, 3]
inline af::dim4 dims_vector(Mesh mesh) { return af::dim4(mesh.nx, mesh.ny, mesh.nz, 3); }
} // namespace magnumafcpp
