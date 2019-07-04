#pragma once
#include "arrayfire.h"
#include <iostream>
struct NonequispacedMesh{
    NonequispacedMesh(){};
    NonequispacedMesh (int nx, int ny, double dx, double dy, std::vector<double> z_spacing);

    //TODO should be const, conficts with string method
    int nx, ny, nz;           //!< Number of cells in x, y, z
    double dx, dy;            //!< Distance between equidistant x, y cells
    std::vector<double> z_spacing;   //
    int nx_expanded, ny_expanded; // Expanded cell sizes for demag FFT
    af::dim4 dims;

    void print(std::ostream& stream);
};
