#pragma once
#include "arrayfire.h"
#include <iostream>
#include <array>
#include <cstdint>

namespace magnumafcpp{

struct Mesh{
    Mesh (uint32_t, uint32_t, uint32_t, double, double, double);
    Mesh (){};

    uint32_t n0, n1, n2;               // Number of cells in x, y, z
    double dx, dy, dz;            // Distance between cells
    double V;                   // Volume of one cell
    int n0_exp, n1_exp, n2_exp; // Expanded cell sizes for demag FFT
    af::dim4 dims;
    af::dim4 dims_expanded;
    void print(std::ostream& stream);
    af::array skyrmconf(const bool point_up = false);
    af::array ellipse(std::array<double, 3> vector, const bool verbose = true);
    af::array ellipse(const int xyz = 0, const bool positive_direction = true);
    af::array init_vortex(const bool positive_direction = true);
    af::array init_sp4();
};
}// namespace magnumafcpp
