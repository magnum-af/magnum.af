#include "mesh.hpp"

namespace magnumafcpp {

Mesh::Mesh(unsigned nx, unsigned ny, unsigned nz, double dx, double dy, double dz)
    : n0(nx), n1(ny), n2(nz), dx(dx), dy(dy), dz(dz), n0_exp(2 * n0), n1_exp(2 * n1), n2_exp((n2 == 1) ? 1 : 2 * n2) {}

std::ostream& operator<<(std::ostream& os, const Mesh& mesh) {
    os << "n0=" << mesh.n0 << " n1=" << mesh.n1 << " n2=" << mesh.n2 << " dx=" << mesh.dx << " dy=" << mesh.dy
       << " dz=" << mesh.dz << " n0_exp=" << mesh.n0_exp << " n1_exp=" << mesh.n1_exp << " n2_exp=" << mesh.n2_exp;
    return os;
}
} // namespace magnumafcpp
