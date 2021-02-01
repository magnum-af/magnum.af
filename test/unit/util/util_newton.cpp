#include "util/newton_iteration.hpp"
#include <cmath>
#include <gtest/gtest.h>

using namespace magnumafcpp;

double x_sqare(double x) { return pow(x, 2) - 1; }

TEST(Util, NetwonIteration) {
    NewtonIteration ni(x_sqare, Verbose(false));
    auto res = ni.run(X0(1e10), Precision(1e-20), EpsilonFactor(1e-4), Imax(100));
    EXPECT_NEAR(res.first, 1, 1e-20);
    EXPECT_NEAR(res.second, 0, 1e-20);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
