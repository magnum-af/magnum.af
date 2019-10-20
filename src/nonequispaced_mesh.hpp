#pragma once
#include "arrayfire.h"
#include <ostream>
#include <iostream>

namespace magnumaf{

struct NonequispacedMesh{
    NonequispacedMesh(){};
    NonequispacedMesh (int nx, int ny, float dx, float dy, std::vector<float> z_spacing);

    //TODO should be const, conficts with string method
    int nx, ny, nz;           //!< Number of cells in x, y, z
    float dx, dy;            //!< Distance between equidistant x, y cells
    std::vector<float> z_spacing;   //
    int nx_expanded, ny_expanded; // Expanded cell sizes for demag FFT
    af::dim4 dims;

    void print(std::ostream& stream = std::cout);
};
}// namespace magnumaf
