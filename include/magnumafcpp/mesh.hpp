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
    af::dim4 dims;                   // Dimension for vector field on mesh, i.e. [nx, ny, nz, 3]
    af::dim4 dims_scalar;            // Dimension for scalar field on mesh, i.e. [nx, ny, nz, 1]
    af::dim4 dims_expanded; // Expanded dimension for demag field, i.e. [2 * nx, 2 * ny, nz == 1 ? 1 : 2 * nz, 3]
    void print(std::ostream& stream);
    af::array skyrmconf(const bool point_up = false);
    af::array ellipse(std::array<double, 3> vector, const bool verbose = true);
    af::array ellipse(const unsigned xyz = 0, const bool positive_direction = true);
    af::array init_vortex(const bool positive_direction = true);
    af::array init_sp4();
};
} // namespace magnumafcpp
