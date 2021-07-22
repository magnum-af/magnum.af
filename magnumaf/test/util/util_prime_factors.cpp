#include "util/prime_factors.hpp"
#include <gmock/gmock.h>
#include <gtest/gtest.h>

using namespace magnumafcpp::util;

TEST(Util, prime_factors) {
    EXPECT_THAT(prime_factors(0), testing::ElementsAre(0));
    EXPECT_THAT(prime_factors(1), testing::ElementsAre(1));
    EXPECT_THAT(prime_factors(2), testing::ElementsAre(2));
    EXPECT_THAT(prime_factors(3), testing::ElementsAre(3));
    EXPECT_THAT(prime_factors(12), testing::ElementsAre(2, 2, 3));
    EXPECT_THAT(prime_factors(40), testing::ElementsAre(2, 2, 2, 5));
    EXPECT_THAT(prime_factors(120), testing::ElementsAre(2, 2, 2, 3, 5));
    EXPECT_THAT(prime_factors(12345), testing::ElementsAre(3, 5, 823));
}

TEST(Util, max_of_prime_factors) {
    EXPECT_EQ(max_of_prime_factors(0), 0);
    EXPECT_EQ(max_of_prime_factors(1), 1);
    EXPECT_EQ(max_of_prime_factors(2), 2);
    EXPECT_EQ(max_of_prime_factors(12), 3);
    EXPECT_EQ(max_of_prime_factors(40), 5);
    EXPECT_EQ(max_of_prime_factors(120), 5);
    EXPECT_EQ(max_of_prime_factors(12345), 823);
}
int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
