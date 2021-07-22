//#include "field_terms/nonequi/nonequi_term.hpp"
//#include "nonequispaced_mesh.hpp"
#include "nonequi/nonequi_demag_field.hpp"
#include "state.hpp"
#include <gtest/gtest.h>

using namespace magnumaf;
TEST(State, integral_nonequimesh_dz_as_gauss_sum) {
    const double dx = 1e-9, dy = 2e-9, dz = 3e-9;
    const int nx = 1, ny = 1, nz = 10;
    std::vector<double> z_spacing;
    for (int i = 1; i <= nz; i++) {
        z_spacing.push_back(i * dz); // gauss sum: dz * (1 + 2 + 3 + ... + n)
    }
    NonequiMesh mesh(nx, ny, dx, dy, z_spacing);
    af::array m = af::constant(0.0, nemesh::dims_v(mesh), f64);
    for (int i = 0; i < nz; i++) {
        if (i % 2 == 0) {
            m(af::span, af::span, i, 2) = 1;
        } else {
            m(af::span, af::span, i, 2) = -1;
        }
    }
    const double Ms = 1.;
    State state(m, Ms);
    // calculating gauss sums s_n = sum_{i=1}^n (a_1 + (i-1)d) = n \frac{a_1 +
    // a_n}{2}
    double positive_sum = nz / 2 * (1 + (nz - 1)) / 2; // 1 + 3 + 5 ... + 9
    double negative_sum = -nz / 2 * (2 + nz) / 2;      // 2 + 4 + 6 ... + 10
    double result = dx * nx * dy * ny * dz * (positive_sum + negative_sum);
    NonequiDemagField nedemag(mesh);
    EXPECT_NEAR((-2. / constants::mu0) * nedemag.Energy_in_J(state, af::constant(1.0, nemesh::dims_v(mesh), f64)), result,
                1e-35);
    // If integral_nonequimesh were public, we could test like this:
    // EXPECT_NEAR((-2. / constants::mu0) * nedemag.integral_nonequimesh(m,
    // state), result, 1e-35);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
