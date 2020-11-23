#include "nonequispaced_mesh.hpp"
#include "func.hpp"

namespace magnumafcpp {

NonequispacedMesh::NonequispacedMesh(unsigned nx, unsigned ny, double dx, double dy, std::vector<double> z_spacing)
    : nx(nx), ny(ny), nz(z_spacing.size()), dx(dx), dy(dy), z_spacing(z_spacing), nx_expanded(2 * nx),
      ny_expanded(2 * ny) {}

std::ostream& operator<<(std::ostream& os, const NonequispacedMesh& nemesh) {
    // void NonequispacedMesh::print(std::ostream& stream) {
    os << "nx=" << nemesh.nx << " ny=" << nemesh.ny << " nz=" << nemesh.nz << " dx=" << nemesh.dx << " dy=" << nemesh.dy
       << " nx_expanded=" << nemesh.nx_expanded << " ny_expanded=" << nemesh.ny_expanded << std::endl;

    os << " dz: ";
    for (auto const& dz : nemesh.z_spacing) {
        os << dz << " ";
    }
    return os;
}
} // namespace magnumafcpp
