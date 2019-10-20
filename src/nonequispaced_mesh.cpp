#include "nonequispaced_mesh.hpp"
#include "func.hpp"

namespace magnumaf{


NonequispacedMesh::NonequispacedMesh (int nx, int ny, float dx, float dy, std::vector<float> z_spacing):
    nx(nx), ny(ny), nz((int) z_spacing.size()), dx(dx), dy(dy), z_spacing(z_spacing), nx_expanded(2*nx), ny_expanded(2*ny), dims(af::dim4(nx, ny, nz, 3))
{
}

void NonequispacedMesh::print(std::ostream& stream){
    stream << "nx=" << nx << " ny=" << ny << " nz=" << nz << " dx=" << dx \
        << " dy=" << dy << " nx_expanded=" << nx_expanded << " ny_expanded=" \
        << ny_expanded << " dims=" << dims << std::endl;

    stream << " dz: ";
    for (auto const& dz : z_spacing){
        stream  << dz << " ";
    }
    stream << std::endl;
}
}// namespace magnumaf
