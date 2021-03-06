#include "util/util.hpp"
#include "util/host_ptr_accessor.hpp"
#include <gtest/gtest.h>

using namespace magnumaf;

TEST(Util, SerialTriangularMatrixTest) {
    // serialized index k, matrix indices i, j
    for (int n = 1; n < 20; n++) {
        std::vector<double> values;
        int k = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i <= j) {
                    values.push_back(0);
                    const int kk = util::ij2k(i, j, n);
                    std::pair<int, int> ij = util::k2ij(k, n);

                    // compare calculated indices with loop indices
                    EXPECT_EQ(i, ij.first);
                    EXPECT_EQ(j, ij.second);
                    EXPECT_EQ(k, kk);

                    // lamda wrapping (just for demonstration)
                    auto k2ij_lambda = [n](const int k) -> std::pair<int, int> { return util::k2ij(k, n); };
                    auto ij2k_lambda = [n](const int i, const int j) -> int { return util::ij2k(i, j, n); };

                    const int kk_l = ij2k_lambda(i, j);
                    std::pair<int, int> ij_l = k2ij_lambda(k);

                    EXPECT_EQ(i, ij_l.first);
                    EXPECT_EQ(j, ij_l.second);
                    EXPECT_EQ(k, kk_l);

                    k++;
                }
            }
        }
        // check number of elements
        EXPECT_EQ(values.size(), (n * (n + 1)) / 2);
    }
}

TEST(Util, 2DimStrideAccessTest) {
    const int ni = 5, nj = 4, nk = 1, nl = 1;
    af::array a = af::iota(af::dim4(ni, nj, nk, nl), af::dim4(1, 1, 1, 1), f64);

    util::HostPtrAccessor<double> host(a);
    for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {
            EXPECT_EQ(util::afvalue_as_f64(a(i, j)), util::stride(i, j, ni));
        }
    }
}

TEST(Util, 3DimStrideAccessTest) {
    const int ni = 5, nj = 4, nk = 3, nl = 1;
    af::array a = af::iota(af::dim4(ni, nj, nk, nl), af::dim4(1, 1, 1, 1), f64);

    util::HostPtrAccessor<double> host(a);
    for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {
            for (int k = 0; k < nk; k++) {
                EXPECT_EQ(util::afvalue_as_f64(a(i, j, k)), util::stride(i, j, k, ni, nj));
            }
        }
    }
}

TEST(Util, 4DimStrideAccessTest) {
    const int ni = 5, nj = 4, nk = 3, nl = 2;
    af::array a = af::iota(af::dim4(ni, nj, nk, nl), af::dim4(1, 1, 1, 1), f64);

    util::HostPtrAccessor<double> host(a);

    for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {
            for (int k = 0; k < nk; k++) {
                for (int l = 0; l < nl; l++) {
                    EXPECT_EQ(util::afvalue_as_f64(a(i, j, k, l)), util::stride(i, j, k, l, ni, nj, nk));
                }
            }
        }
    }
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
