#include "field_terms/nonequi/nonequi_external_field.hpp"
#include "field_terms/micro/external_field.hpp"
#include "nonequispaced_mesh.hpp"
#include <gtest/gtest.h>
#include <vector>

using namespace magnumaf;

// TEST(NonequiExternalField, E) {
//    int nx = 1;
//    int ny = 2;
//    double dx = 1;
//    double dy = 2;
//    std::vector<double> dz = {1, 2, 3};
//    NonequiMesh nemesh(nx, ny, dx, dy, dz);
//    af::array m = af::constant(1., nx, ny, nemesh.nz, 3, f64);
//    af::array field = af::constant(1., nx, ny, nemesh.nz, 3, f64);
//    NonequiExternalField nonequi_external(nemesh, field);
//    ExternalField external(field);
//    double Ms = 1;
//    State state(m, Ms);
//    // TODO:
//    // EXPECT_EQ(nonequi_external.Energy_in_J(state_nemesh), external.Energy_in_J(state_nemesh));
//    // State state_mesh(Mesh(nx, ny, nemesh.nz, dx, dy, dz[0]), 1, m);
//    // EXPECT_EQ(nonequi_external.Energy_in_J(state_nemesh), external.Energy_in_J(state_mesh));
//    // EXPECT_EQ(nonequi_external.Energy_in_J(state_nemesh), 0);
//}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
