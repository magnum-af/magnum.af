#include "nonequispaced_mesh.hpp"
#include "func.hpp"

namespace magnumafcpp {

NonequiMesh::NonequiMesh(unsigned nx, unsigned ny, double dx, double dy, std::vector<double> z_spacing)
    : nx(nx), ny(ny), nz(z_spacing.size()), dx(dx), dy(dy), z_spacing(z_spacing) {}

std::ostream& operator<<(std::ostream& os, const NonequiMesh& nemesh) {
    os << "nx=" << nemesh.nx << " ny=" << nemesh.ny << " nz=" << nemesh.nz << " dx=" << nemesh.dx << " dy=" << nemesh.dy
       << " dz: ";
    for (auto const& dz : nemesh.z_spacing) {
        os << dz << " ";
    }
    return os;
}

} // namespace magnumafcpp
