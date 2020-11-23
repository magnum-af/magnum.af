#include "mesh.hpp"

namespace magnumafcpp {

Mesh::Mesh(unsigned nx, unsigned ny, unsigned nz, double dx, double dy, double dz)
    : n0(nx), n1(ny), n2(nz), dx(dx), dy(dy), dz(dz) {}

std::ostream& operator<<(std::ostream& os, const Mesh& mesh) {
    os << "n0=" << mesh.n0 << " n1=" << mesh.n1 << " n2=" << mesh.n2 << " dx=" << mesh.dx << " dy=" << mesh.dy
       << " dz=" << mesh.dz;
    return os;
}
} // namespace magnumafcpp
