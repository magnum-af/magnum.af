#include "state.hpp"
#include "func.hpp"
#include <gtest/gtest.h>
#include <vector>

using namespace magnumafcpp;
TEST(State, Init_Ms_field) {
    int nx = 1, ny = 1, nz = 1;
    auto m = af::constant(0, nx, ny, nz, 3, f64);
    m(af::span, af::span, af::span, 0) = 1;
    auto Ms_scalar = af::constant(1, nx, ny, nz, 1, f64);
    State state_scalar(Mesh(1, 1, 1, 1, 1, 1), Ms_scalar, m);
    EXPECT_EQ(state_scalar.Ms_field.dims(3), 3);

    auto Ms_vector = af::constant(1, nx, ny, nz, 3, f64);
    State state_legacy(Mesh(1, 1, 1, 1, 1, 1), Ms_vector, m);
    EXPECT_EQ(state_legacy.Ms_field.dims(3), 3);

    // To test wrapping we use python, elso somting like this:
    // auto Ms_wrapping = af::constant(1, nx, ny, nz, 1, f64);
    // State state_wrapping(Mesh(1,1,1,1,1,1), (long int) Ms_scalar.get(), (long
    // int) m.get()); EXPECT_EQ(state_wrapping.Ms_field.dims(3), 3);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
