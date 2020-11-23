#pragma once
#include "arrayfire.h"
#include <iostream>

namespace magnumafcpp {

struct Mesh {
    Mesh(unsigned nx, unsigned ny, unsigned nz, double dx, double dy, double dz)
        : n0(nx), n1(ny), n2(nz), dx(dx), dy(dy), dz(dz) {}

    unsigned n0, n1, n2;             // Number of cells in x, y, z
    double dx, dy, dz;               // Distance between cells
};

inline std::ostream& operator<<(std::ostream& os, const Mesh& mesh) {
    os << "n0=" << mesh.n0 << " n1=" << mesh.n1 << " n2=" << mesh.n2 << " dx=" << mesh.dx << " dy=" << mesh.dy
       << " dz=" << mesh.dz;
    return os;
}

// Dimension for scalar field on mesh, i.e. [nx, ny, nz, 1]
inline af::dim4 dims_scalar(Mesh mesh) { return af::dim4(mesh.n0, mesh.n1, mesh.n2, 1); }

// Dimension for vector field on mesh, i.e. [nx, ny, nz, 3]
inline af::dim4 dims_vector(Mesh mesh) { return af::dim4(mesh.n0, mesh.n1, mesh.n2, 3); }
} // namespace magnumafcpp
