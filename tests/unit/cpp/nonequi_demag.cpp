#include <gtest/gtest.h>
#include "../../../src/mesh.cpp"
#include "../../../src/func.cpp"
#include "../../../src/misc.cpp"
#include "../../../src/llg_terms/micro_demag_nonequi.cpp"
#include "../../../src/llg_terms/micro_demag.cpp"
 
// Exemplary unit test
TEST(NonEquiDemagNxxNxyNearTest, n) {
    double Nxx_non = newell_nonequi::Nxx_nonequi(1, 2, 3, 1, 2, 3, 1, 2, 3);
    double Nxx_f = newell::Nxx(1, 2, 3, 1, 2, 3);
    double Nxy_non = newell_nonequi::Nxy_nonequi(1, 2, 3, 1, 2, 3, 1, 2, 3);
    double Nxyg = newell::Nxy(1, 2, 3, 1, 2, 3);

    EXPECT_NEAR(Nxx_non, Nxx_f, 1e-14);
    double Nxx_rel_diff = fabs(2*(Nxx_non-Nxx_f)/(Nxx_non+Nxx_f));
    EXPECT_NEAR(Nxx_rel_diff, 0, 1e-11);

    EXPECT_NEAR(Nxy_non, Nxyg, 1e-15);
    double Nxy_rel_diff = fabs(2*(Nxy_non-Nxyg)/(Nxy_non+Nxyg));
    EXPECT_NEAR(Nxy_rel_diff, 0, 1e-10);
}

TEST(NonEquiDemagNxxNxyFarTest, n) {
    double Nxx_non = newell_nonequi::Nxx_nonequi(1, 2, 300, 1, 2, 3, 1, 2, 3);
    double Nxx_f = newell::Nxx(1, 2, 300, 1, 2, 3);
    double Nxy_non = newell_nonequi::Nxy_nonequi(1, 2, 300, 1, 2, 3, 1, 2, 3);
    double Nxyg = newell::Nxy(1, 2, 300, 1, 2, 3);

    double Nxx_rel_diff = fabs(2*(Nxx_non-Nxx_f)/(Nxx_non+Nxx_f));
    double Nxx_abs_diff = fabs(Nxx_non-Nxx_f);

    EXPECT_NEAR(Nxx_rel_diff, 0, 1e1);
    EXPECT_NEAR(Nxx_abs_diff, 0, 1e-8);

    double Nxy_rel_diff = fabs(2*(Nxy_non-Nxyg)/(Nxy_non+Nxyg));
    double Nxy_abs_diff = fabs(Nxy_non-Nxyg);

    EXPECT_NEAR(Nxy_rel_diff, 0, 1e1);
    EXPECT_NEAR(Nxy_abs_diff, 0, 1e-12);
}

TEST(NonEquiDemagUnitCubeTest, n) {
    const double a = 1;
    const int  ix = 2;
    double Nxx = newell::Nxx(ix, 0, 0, a, a, a);
    double Nxx_ne = newell_nonequi::Nxx_nonequi(ix, 0, 0, a, a, a, a, a, a);
    EXPECT_NEAR(Nxx, Nxx_ne, 1e-16);
}

//Testing layout in x: #|##
TEST(NonEquiDemagTwoCubesVersusOneCubeTest_x, n) {
    const double a = 1;
    const int ix = 1;
    const double Nxx_1 = newell::Nxx(  ix, 0, 0, a, a, a);
    const double Nxx_2 = newell::Nxx(2*ix, 0, 0, a, a, a);
    const double Nxx = Nxx_1 + Nxx_2;
    const double Nxx_ne = newell_nonequi::Nxx_nonequi(ix, 0, 0, a, a, a, 2*a, a, a);
    EXPECT_NEAR(Nxx, Nxx_ne, 1e-16);
}

//Testing layout in y: #|##
TEST(NonEquiDemagTwoCubesVersusOneCubeTest_y, n) {
    const double a = 1;
    const int iy = 1;
    const double Nxx_1 = newell::Nxx(0,   iy, 0, a, a, a);
    const double Nxx_2 = newell::Nxx(0, 2*iy, 0, a, a, a);
    const double Nxx = Nxx_1 + Nxx_2;
    const double Nxx_ne = newell_nonequi::Nxx_nonequi(0, iy, 0, a, a, a, a, 2*a, a);
    EXPECT_NEAR(Nxx, Nxx_ne, 1e-15);
}

//Testing layout in z: #|##
TEST(NonEquiDemagTwoCubesVersusOneCubeTest_z, n) {
    const double a = 1;
    const int iz = 1;
    const double Nxx_1 = newell::Nxx(0,   iz, 0, a, a, a);
    const double Nxx_2 = newell::Nxx(0, 2*iz, 0, a, a, a);
    const double Nxx = Nxx_1 + Nxx_2;
    const double Nxx_ne = newell_nonequi::Nxx_nonequi(0, iz, 0, a, a, a, a, 2*a, a);
    EXPECT_NEAR(Nxx, Nxx_ne, 1e-15);
}

//                      |#
//Testing layout in x: #|#
TEST(NonEquiDemagTwoCubesVersusOneCubeTest_xy, n) {
    const double a = 1;
    const int ixy = 1;
    const double Nxx_1 = newell::Nxx(ixy,   ixy, 0, a, a, a);
    const double Nxx_2 = newell::Nxx(ixy, 2*ixy, 0, a, a, a);
    const double Nxx = Nxx_1 + Nxx_2;
    const double Nxx_ne = newell_nonequi::Nxx_nonequi(ixy, ixy, 0, a, a, a, a, 2*a, a);
    EXPECT_NEAR(Nxx, Nxx_ne, 1e-15);
}

//Testing layout in x: #|#
//                      |#
TEST(NonEquiDemagTwoCubesVersusOneCubeTest_x_minus_y, n) {
    const double a = 1;
    const int ixy = 1;
    const double Nxx_1 = newell::Nxx(ixy, -  ixy, 0, a, a, a);
    const double Nxx_2 = newell::Nxx(ixy, -2*ixy, 0, a, a, a);
    const double Nxx = Nxx_1 + Nxx_2;
    const double Nxx_ne = newell_nonequi::Nxx_nonequi(ixy, -2*ixy, 0, a, a, a, a, 2*a, a);
    EXPECT_NEAR(Nxx, Nxx_ne, 1e-15);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
