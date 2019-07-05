#include <gtest/gtest.h>
#include "../../../src/util/util.hpp"
#include <cmath>

using namespace magnumaf;


// calculating arithmetic series a_n = n i.e:  {1, 2, ..., n-1, n}
// standard deviation following http://clay6.com/qa/49501/determine-mean-and-standard-deviation-of-first-n-terms-of-an-a-p-whose-firs
TEST(Util, mean_stdev_no_minus) {
    const int nn = 10000;
    std::vector<double> vec;
    for(int i = 1; i<=nn; i++){
        vec.push_back(i);
    }
    auto res = mean_stdev_no_minus(vec);
    double n = nn;
    EXPECT_EQ(res.first, n * (n + 1.)/2./n);
    EXPECT_NEAR(res.second, std::sqrt( (std::pow(n, 2) -1) / 12.), 0.);
    EXPECT_NEAR(res.second, std::sqrt( (n+1) * (n-1) / 12.), 0.);// equivalent to prev
}


TEST(Util, mean_stdev_w_minus) {
    const int nn = 10000;
    std::vector<double> vec;
    for(int i = 1; i<=nn; i++){
        vec.push_back(i);
    }
    auto res = mean_stdev_w_minus(vec);
    double n = nn;
    EXPECT_EQ(res.first, n * (n + 1.)/2./n);
    EXPECT_NEAR(res.second, std::sqrt( (std::pow(n, 2) -1) / 12.), 0.15);
    EXPECT_NEAR(res.second, std::sqrt( (n+1) * (n-1) / 12.), 0.15);// equivalent to prev
}


int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
