#include "../../src/mesh.hpp"
#include "../../src/nonequispaced_mesh.hpp"//TODO needed for vtk_IO, include headers instead?
#include "../../src/func.hpp"
#include "../../src/state.hpp"
#include "../../src/misc.hpp"
#include "../../src/vtk_IO.hpp"
#include "../../src/llg_terms/micro_demag_nonequi.hpp"
//#include "../../src/llg_terms/micro_demag_nonequi.cpp"
#include "../../src/llg_terms/micro_demag.hpp"
//#include "../../src/llg_terms/micro_demag.cpp"
#include <gtest/gtest.h>
#include <vector>

using namespace magnumafcpp;


// Exemplary unit test
TEST(NonEquiDemag, NxxNxyNearTest) {
    int ix = 1;
    int iy = 2;
    int iz = 3;
    double dx = 1;
    double dy = 2;
    double dz = 3;
    double dX = dx;
    double dY = dy;
    double dZ = dz;
    double x = (double) ix * dx;
    double y = (double) iy * dy;
    double z = (double) iz * dz;

    double Nxx_non = newell_nonequi::Nxx(x, y, z, dx, dy, dz, dX, dY, dZ)/(dX * dY * dZ);//multiply by 1/tau as tau is now added in Heff calculation
    double Nxx_f = newell::Nxx(ix, iy, iz, dx, dy, dz);
    double Nxy_non = newell_nonequi::Nxy(x, y, z, dx, dy, dz, dX, dY, dZ)/(dX * dY * dZ);
    double Nxyg = newell::Nxy(ix, iy, iz, dx, dy, dz);

    EXPECT_NEAR(Nxx_non, Nxx_f, 1e-14);
    double Nxx_rel_diff = fabs(2*(Nxx_non-Nxx_f)/(Nxx_non+Nxx_f));
    EXPECT_NEAR(Nxx_rel_diff, 0, 1e-11);

    EXPECT_NEAR(Nxy_non, Nxyg, 1e-15);
    double Nxy_rel_diff = fabs(2*(Nxy_non-Nxyg)/(Nxy_non+Nxyg));
    EXPECT_NEAR(Nxy_rel_diff, 0, 1e-10);
}

TEST(NonEquiDemag, NxxNxyFarTest) {
    int ix = 1;
    int iy = 2;
    int iz = 300;
    double dx = 1;
    double dy = 2;
    double dz = 3;
    double dX = dx;
    double dY = dy;
    double dZ = dz;
    double x = (double) ix * dx;
    double y = (double) iy * dy;
    double z = (double) iz * dz;

    double Nxx_non = newell_nonequi::Nxx(x, y, z, dx, dy, dz, dX, dY, dZ)/(dX * dY * dZ);
    double Nxx_f = newell::Nxx(ix, iy, iz, dx, dy, dz);
    double Nxy_non = newell_nonequi::Nxy(x, y, z, dx, dy, dz, dX, dY, dZ)/(dX * dY * dZ);
    double Nxyg = newell::Nxy(ix, iy, iz, dx, dy, dz);

    double Nxx_rel_diff = fabs(2*(Nxx_non-Nxx_f)/(Nxx_non+Nxx_f));
    double Nxx_abs_diff = fabs(Nxx_non-Nxx_f);

    EXPECT_NEAR(Nxx_rel_diff, 0, 1e1);
    EXPECT_NEAR(Nxx_abs_diff, 0, 1e-8);

    double Nxy_rel_diff = fabs(2*(Nxy_non-Nxyg)/(Nxy_non+Nxyg));
    double Nxy_abs_diff = fabs(Nxy_non-Nxyg);

    EXPECT_NEAR(Nxy_rel_diff, 0, 1e1);
    EXPECT_NEAR(Nxy_abs_diff, 0, 1e-12);
}

TEST(NonEquiDemag, DistanceFromIndexTest) {
    std::vector<double> vz = {1, 2, 3};
    EXPECT_EQ(newell_nonequi::nonequi_index_distance(vz, 0, 1), 1);
    EXPECT_EQ(newell_nonequi::nonequi_index_distance(vz, 0, 2), 1 + 2);
    EXPECT_EQ(newell_nonequi::nonequi_index_distance(vz, 1, 0), -1);
    EXPECT_EQ(newell_nonequi::nonequi_index_distance(vz, 2, 0), -1 - 2);
    EXPECT_EQ(newell_nonequi::nonequi_index_distance(vz, 0, 3, false), 1 + 2 + 3);//Note:
    //This is not out of bound because distance is from i to j-1. The arising warning is surpressed here.
}

TEST(NonEquiDemag, UnitCubeTest) {
    const double a = 2;
    const int  ix = 2;
    double Nxx = newell::Nxx(ix, 0, 0, a, a, a);
    double Nxx_ne = newell_nonequi::Nxx(ix * a, 0, 0, a, a, a, a, a, a)/(a * a * a);
    EXPECT_NEAR(Nxx, Nxx_ne, 1e-16);

    double Nxy = newell::Nxy(ix, 0, 0, a, a, a);
    double Nxy_ne = newell_nonequi::Nxy(ix * a, 0, 0, a, a, a, a, a, a)/(a * a * a);
    EXPECT_EQ(Nxy, 0);
    EXPECT_EQ(Nxy_ne, 0);
}

TEST(NonEquiDemag, SymmetryTest) {
    const double dx = 2, dy = 3, dz = 4;
    const double dX = 3, dY = 4, dZ = 5;
    const double x = 20, y = 22, z = 24;

    double Nxx = newell_nonequi::Nxx(x, y, z, dx, dy, dz, dX, dY, dZ);
    double Nxx_switch = newell_nonequi::Nxx(-x, -y, -z, dX, dY, dZ, dx, dy, dz);
    EXPECT_NEAR(Nxx, Nxx_switch, 1e-11);
}

//Testing layout in x: #|##
TEST(NonEquiDemag, TwoCubesVersusOneCubeTest_x) {
    const double a = 2;
    const int ix = 1;

    const double Nxx_1 = newell::Nxx(  ix, 0, 0, a, a, a);
    const double Nxx_2 = newell::Nxx(2*ix, 0, 0, a, a, a);
    const double Nxx = Nxx_1 + Nxx_2;
    const double Nxx_ne = newell_nonequi::Nxx(- ix * a, 0, 0, 2*a, a, a, a, a, a)/(a * a * a);
    const double Nxx_ne2 = newell_nonequi::Nxx(2 * ix * a, 0, 0, 2*a, a, a, a, a, a)/(a * a * a);// x-symmetric case
    EXPECT_NEAR(Nxx, Nxx_ne, 1e-16);
    EXPECT_NEAR(Nxx, Nxx_ne2, 1e-16);

    const double Nxy_1 = newell::Nxy(  ix, 0, 0, a, a, a);
    const double Nxy_2 = newell::Nxy(2*ix, 0, 0, a, a, a);
    const double Nxy = Nxy_1 + Nxy_2;
    const double Nxy_ne = newell_nonequi::Nxy(- ix * a, 0, 0, 2*a, a, a, a, a, a)/(a * a * a);
    const double Nxy_ne2 = newell_nonequi::Nxy(2 * ix * a, 0, 0, 2*a, a, a, a, a, a)/(a * a * a);// x-symmetric case
    EXPECT_NEAR(Nxy, Nxy_ne, 1e-16);
    EXPECT_NEAR(Nxy, Nxy_ne2, 1e-16);
}

//Testing layout in y: #|##
TEST(NonEquiDemag, TwoCubesVersusOneCubeTest_y) {
    const double a = 2;
    const int iy = 1;

    const double Nxx_1 = newell::Nxx(0,   iy, 0, a, a, a);
    const double Nxx_2 = newell::Nxx(0, 2*iy, 0, a, a, a);
    const double Nxx = Nxx_1 + Nxx_2;
    const double Nxx_ne = newell_nonequi::Nxx(0, -iy * a, 0, a, 2*a, a, a, a, a)/(a * a * a);
    EXPECT_NEAR(Nxx, Nxx_ne, 1e-15);

    const double Nxy_1 = newell::Nxy(0,   iy, 0, a, a, a);
    const double Nxy_2 = newell::Nxy(0, 2*iy, 0, a, a, a);
    const double Nxy = Nxy_1 + Nxy_2;
    const double Nxy_ne = newell_nonequi::Nxy(0, -iy * a, 0, a, 2*a, a, a, a, a)/(a * a * a);
    EXPECT_NEAR(Nxy, Nxy_ne, 1e-15);
}

//Testing layout in z: #|##
TEST(NonEquiDemag, TwoCubesVersusOneCubeTest_z) {
    const double a = 2;
    const int iz = 1;

    const double Nxx_1 = newell::Nxx(0, 0,   iz, a, a, a);
    const double Nxx_2 = newell::Nxx(0, 0, 2*iz, a, a, a);
    const double Nxx = Nxx_1 + Nxx_2;
    const double Nxx_ne = newell_nonequi::Nxx(0, 0, -iz * a, a, a, 2*a, a, a, a)/(a * a * a);
    EXPECT_NEAR(Nxx, Nxx_ne, 1e-15);
}

//                             |#
//Testing layout in x and +y: #|#
TEST(NonEquiDemag, TwoCubesVersusOneCubeTest_xy) {
    const double a = 2;
    const int ixy = 1;

    const double Nxx_1 = newell::Nxx(ixy,   ixy, 0, a, a, a);
    const double Nxx_2 = newell::Nxx(ixy, 2*ixy, 0, a, a, a);
    const double Nxx = Nxx_1 + Nxx_2;
    const double Nxx_ne = newell_nonequi::Nxx(-ixy * a, -ixy * a, 0, a, 2*a, a, a, a, a)/(a * a * a);
    EXPECT_NEAR(Nxx, Nxx_ne, 1e-15);

    const double Nxy_1 = newell::Nxy(ixy,   ixy, 0, a, a, a);
    const double Nxy_2 = newell::Nxy(ixy, 2*ixy, 0, a, a, a);
    const double Nxy = Nxy_1 + Nxy_2;
    const double Nxy_ne = newell_nonequi::Nxy(-ixy * a, -ixy * a, 0, a, 2*a, a, a, a, a)/(a * a * a);
    EXPECT_NEAR(Nxy, Nxy_ne, 1e-15);
}

//Testing layout in x and -y: #|#
//                             |#
TEST(NonEquiDemag, TwoCubesVersusOneCubeTest_x_minus_y) {
    const double a = 2;
    const int ixy = 1;

    const double Nxx_1 = newell::Nxx(ixy, -  ixy, 0, a, a, a);
    const double Nxx_2 = newell::Nxx(ixy, -2*ixy, 0, a, a, a);
    const double Nxx = Nxx_1 + Nxx_2;
    const double Nxx_ne = newell_nonequi::Nxx(-ixy * a, 2*ixy * a, 0, a, 2*a, a, a, a, a)/(a * a * a);
    EXPECT_NEAR(Nxx, Nxx_ne, 1e-15);

    const double Nxy_1 = newell::Nxy(ixy, -  ixy, 0, a, a, a);
    const double Nxy_2 = newell::Nxy(ixy, -2*ixy, 0, a, a, a);
    const double Nxy = Nxy_1 + Nxy_2;
    const double Nxy_ne = newell_nonequi::Nxy(-ixy * a, 2*ixy * a, 0, a, 2*a, a, a, a, a)/(a * a * a);
    EXPECT_NEAR(Nxy, Nxy_ne, 1e-15);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
