#include <gtest/gtest.h>
#include "../../../src/mesh.cpp"
#include "../../../src/func.cpp"
#include "../../../src/misc.cpp"
#include "../../../src/llg_terms/micro_nonequi_demag.cpp"
 
// Exemplary unit test
TEST(NonEquiDemagNxxNxyNearTest, n) {
    double Nxx_non = newell::Nxx_nonequi(1, 2, 3, 1, 2, 3, 1, 2, 3);
    double Nxxf = newell::Nxxf(1, 2, 3, 1, 2, 3);
    double Nxy_non = newell::Nxy_nonequi(1, 2, 3, 1, 2, 3, 1, 2, 3);
    double Nxyg = newell::Nxxg(1, 2, 3, 1, 2, 3);

    double Nxx_rel_diff = fabs(2*(Nxx_non-Nxxf)/(Nxx_non+Nxxf));
    double Nxx_abs_diff = fabs(Nxx_non-Nxxf);

    EXPECT_NEAR(Nxx_rel_diff, 0, 1e-11);
    EXPECT_NEAR(Nxx_abs_diff, 0, 1e-14);

    double Nxy_rel_diff = fabs(2*(Nxy_non-Nxyg)/(Nxy_non+Nxyg));
    double Nxy_abs_diff = fabs(Nxy_non-Nxyg);

    EXPECT_NEAR(Nxy_rel_diff, 0, 1e-10);
    EXPECT_NEAR(Nxy_abs_diff, 0, 1e-15);
}

TEST(NonEquiDemagNxxNxyFarTest, n) {
    double Nxx_non = newell::Nxx_nonequi(1, 2, 300, 1, 2, 3, 1, 2, 3);
    double Nxxf = newell::Nxxf(1, 2, 300, 1, 2, 3);
    double Nxy_non = newell::Nxy_nonequi(1, 2, 300, 1, 2, 3, 1, 2, 3);
    double Nxyg = newell::Nxxg(1, 2, 300, 1, 2, 3);

    double Nxx_rel_diff = fabs(2*(Nxx_non-Nxxf)/(Nxx_non+Nxxf));
    double Nxx_abs_diff = fabs(Nxx_non-Nxxf);

    EXPECT_NEAR(Nxx_rel_diff, 0, 1e1);
    EXPECT_NEAR(Nxx_abs_diff, 0, 1e-8);

    double Nxy_rel_diff = fabs(2*(Nxy_non-Nxyg)/(Nxy_non+Nxyg));
    double Nxy_abs_diff = fabs(Nxy_non-Nxyg);

    EXPECT_NEAR(Nxy_rel_diff, 0, 1e1);
    EXPECT_NEAR(Nxy_abs_diff, 0, 1e-12);

    //std::cout << "Nxy_nonequi: " << Nxy_non << std::endl;
    //std::cout << "Nxy        : " << Nxyg << std::endl;
    //std::cout << "Nxy_diff   : " << Nxx_rel_diff << std::endl << std::endl;
    //std::cout << "Nxy_diff   : " << Nxx_abs_diff << std::endl << std::endl;
    //std::cout << "Nxy_diff   : " << Nxy_rel_diff << std::endl << std::endl;
    //std::cout << "Nxy_diff   : " << Nxy_abs_diff << std::endl << std::endl;
    //std::cout << "Nxx_nonequi: " << Nxx_non << std::endl;
    //std::cout << "Nxx        : " << Nxxf << std::endl;
    //std::cout << "Nxx_diff   : " << 2*(Nxx_non-Nxxf)/(Nxx_non+Nxxf) << std::endl << std::endl;
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
