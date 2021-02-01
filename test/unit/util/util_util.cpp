#include "util/util.hpp"
#include <cmath>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

using namespace magnumafcpp;

// calculating arithmetic series a_n = n i.e:  {1, 2, ..., n-1, n}
// standard deviation following
// http://clay6.com/qa/49501/determine-mean-and-standard-deviation-of-first-n-terms-of-an-a-p-whose-firs
TEST(Util, mean_stdev_no_minus) {
    const int nn = 10000;
    std::vector<double> vec;
    for (int i = 1; i <= nn; i++) {
        vec.push_back(i);
    }
    auto res = mean_stdev_no_minus(vec);
    double n = nn;
    EXPECT_EQ(res.first, n * (n + 1.) / 2. / n);
    EXPECT_NEAR(res.second, std::sqrt((std::pow(n, 2) - 1) / 12.), 0.);
    EXPECT_NEAR(res.second, std::sqrt((n + 1) * (n - 1) / 12.),
                0.); // equivalent to prev
}

TEST(Util, mean_stdev_w_minus) {
    const int nn = 10000;
    std::vector<double> vec;
    for (int i = 1; i <= nn; i++) {
        vec.push_back(i);
    }
    auto res = mean_stdev_w_minus(vec);
    double n = nn;
    EXPECT_EQ(res.first, n * (n + 1.) / 2. / n);
    EXPECT_NEAR(res.second, std::sqrt((std::pow(n, 2) - 1) / 12.), 0.15);
    EXPECT_NEAR(res.second, std::sqrt((n + 1) * (n - 1) / 12.),
                0.15); // equivalent to prev
}

TEST(Util, cross_product) {
    EXPECT_THAT(cross_product({1, 0, 0}, {1, 0, 0}), testing::ElementsAre(0, 0, 0));
    EXPECT_THAT(cross_product({1, 0, 0}, {0, 1, 0}), testing::ElementsAre(0, 0, 1));
    EXPECT_THAT(cross_product({2, 0, 0}, {0, 2, 0}), testing::ElementsAre(0, 0, 4));
}

TEST(Util, dot_product) {
    EXPECT_EQ(dot_product({0, 0, 0}, {0, 0, 0}), 0);
    EXPECT_EQ(dot_product({1, 0, 0}, {1, 0, 0}), 1);
    EXPECT_EQ(dot_product({2, 0, 0}, {2, 0, 0}), 4);
    EXPECT_EQ(dot_product({1, 0, 0}, {0, 1, 0}), 0);
    EXPECT_EQ(dot_product({1, 2, 3}, {1, 2, 3}), 14);
    EXPECT_EQ(dot_product({1, 1, 1}, {1, 1, 1}), 3);
}

TEST(Util, vector_norm) {
    EXPECT_EQ(vector_norm({0, 0, 0}), 0);
    EXPECT_EQ(vector_norm({1, 0, 0}), 1);
    EXPECT_EQ(vector_norm({2, 0, 0}), 2);
    EXPECT_EQ(vector_norm({0, 3, 0}), 3);
    EXPECT_EQ(vector_norm({0, 0, 4}), 4);
    EXPECT_EQ(vector_norm({1, 1, 1}), std::sqrt(3));
}

TEST(Util, normalize_vector) {
    EXPECT_THAT(normalize_vector({2, 0, 0}), testing::ElementsAre(1, 0, 0));
    EXPECT_THAT(normalize_vector({0, 2, 0}), testing::ElementsAre(0, 1, 0));
    EXPECT_THAT(normalize_vector({0, 0, 2}), testing::ElementsAre(0, 0, 1));
    EXPECT_THAT(normalize_vector({2, 2, 0}), testing::ElementsAre(1 / std::sqrt(2), 1 / std::sqrt(2), 0));
    EXPECT_THAT(normalize_vector({2, 2, 2}),
                testing::ElementsAre(1 / std::sqrt(3), 1 / std::sqrt(3), 1 / std::sqrt(3)));
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
