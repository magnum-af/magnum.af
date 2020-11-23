#pragma once
#include "arrayfire.h"
#include <array>
#include <cstdint>
#include <iostream>

namespace magnumafcpp {

struct Mesh {
    Mesh(unsigned nx, unsigned ny, unsigned nz, double dx, double dy, double dz);

    unsigned n0, n1, n2;             // Number of cells in x, y, z
    double dx, dy, dz;               // Distance between cells
    unsigned n0_exp, n1_exp, n2_exp; // Expanded cell sizes for demag FFT
};

std::ostream& operator<<(std::ostream& os, const Mesh& mesh);

// Dimension for scalar field on mesh, i.e. [nx, ny, nz, 1]
inline af::dim4 dims_scalar(Mesh mesh) { return af::dim4(mesh.n0, mesh.n1, mesh.n2, 1); }

// Dimension for vector field on mesh, i.e. [nx, ny, nz, 3]
inline af::dim4 dims_vector(Mesh mesh) { return af::dim4(mesh.n0, mesh.n1, mesh.n2, 3); }

} // namespace magnumafcpp
