#include "../../src/util/gradient_decent.hpp"
#include <gtest/gtest.h>
#include <cmath>

using namespace magnumafcpp;

double polyn(double x)
{
    return pow(x, 4) - 3 * pow(x, 3) + 2;
}

TEST(Util, GradientDecent)
{
    GradientDecent gd(polyn, X_start_val(6), Precision(1e-6), Gamma(1e-2), maxIters(1000), Epsilon(1e-9), Verbose(false));
    auto res = gd.minimize();
    EXPECT_NEAR(res.first, 2.25, 4 * 1e-6);
    EXPECT_NEAR(res.second, -1675. / 256., 1e-6);
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
