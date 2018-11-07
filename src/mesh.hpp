#ifndef MESH_H
#define MESH_H
#include <iostream>
#include "arrayfire.h"
struct Mesh{
    int n0,n1,n2;               // Number of cells in x,y,z
    double dx,dy,dz;            // Distance between cells
    double V;                   // Volume of one cell
    int n0_exp, n1_exp, n2_exp; // Expanded cell sizes for demag FFT
    Mesh (int, int, int, double, double, double);
    af::dim4 dims;
    af::dim4 dims_expanded;
    void print(std::ostream& stream);
    af::array skyrmconf(const bool point_up = false);
    af::array ellipse(long int& n_cells, const int xyz = 0, const bool positive_direction = true);
    af::array init_vortex(long int& n_cells, const bool positive_direction = true);
};
#endif
