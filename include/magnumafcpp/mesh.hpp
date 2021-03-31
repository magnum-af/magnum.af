#pragma once
#include "arrayfire.h"
#include <cstddef>
#include <iostream>

//! C++ magnum.af library
namespace magnumafcpp {

using std::size_t;

struct Mesh {
    constexpr Mesh(size_t nx, size_t ny, size_t nz, double dx, double dy, double dz); // get rid of ctor?
    size_t nx, ny, nz;   ///< Number of cells in x, y, z
    double dx, dy, dz;   ///< Length of a cell in x, y, z in [m]
};

// avoiding shadow warning:
inline constexpr Mesh::Mesh(size_t nx_, size_t ny_, size_t nz_, double dx_, double dy_, double dz_)
    : nx(nx_), ny(ny_), nz(nz_), dx(dx_), dy(dy_), dz(dz_) {}

inline std::ostream& operator<<(std::ostream& os, const Mesh& mesh) {
    os << "nx=" << mesh.nx << " ny=" << mesh.ny << " nz=" << mesh.nz << " dx=" << mesh.dx << " dy=" << mesh.dy
       << " dz=" << mesh.dz;
    return os;
}

//! Mesh class helper functions
namespace mesh {
// Dimension for scalar field on mesh, i.e. af::dim4(nx, ny, nz, 1)
inline af::dim4 dims_s(Mesh mesh) { return af::dim4(mesh.nx, mesh.ny, mesh.nz, 1); }

// Dimension for vector field on mesh, i.e. af::dim4(nx, ny, nz, 3)
inline af::dim4 dims_v(Mesh mesh) { return af::dim4(mesh.nx, mesh.ny, mesh.nz, 3); }
} // namespace mesh
} // namespace magnumafcpp
