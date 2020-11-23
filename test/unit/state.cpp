#include "state.hpp"
#include "func.hpp"
#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
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

TEST(State, normalize_m) {
    Mesh mesh(1, 1, 1, 1e-9, 1e-9, 1e-9);
    double Ms = 1;
    // af::array m = af::constant(std::sqrt(3), 1, 1, 1, 3, f64); // TODO check why this throws warning, numerical?
    af::array m = af::constant(1, 1, 1, 1, 3, f64);
    State state(mesh, Ms, m);
    EXPECT_EQ(state.m(0, 0, 0, 0).scalar<double>(), 1 / std::sqrt(3));
    EXPECT_EQ(state.m(0, 0, 0, 1).scalar<double>(), 1 / std::sqrt(3));
    EXPECT_EQ(state.m(0, 0, 0, 2).scalar<double>(), 1 / std::sqrt(3));
}

// Testing State::mean_m() for scalar and array (zeros/nozeros) Ms
TEST(State, mean_m) {
    Mesh mesh(2, 2, 2, 1e-9, 1e-9, 1e-9);
    double Ms = 1e6;

    auto test = [](State state, std::array<double, 3> res) { EXPECT_THAT(state.mean_m(), res); };

    // vec_as_array method:
    auto m = [dims = dims_scalar(mesh)](std::array<double, 3> vec) {
        return af::tile(af::array(af::dim4(1, 1, 1, 3), vec.data()), dims);
    };

    // scalar Ms
    test({mesh, Ms, m({1, 0, 0})}, {1, 0, 0});
    test({mesh, Ms, m({0, 1, 0})}, {0, 1, 0});
    test({mesh, Ms, m({0, 0, 1})}, {0, 0, 1});

    // Ms.field
    auto Ms_field = af::constant(Ms, dims_scalar(mesh), f64);
    test({mesh, Ms_field, m({1, 0, 0})}, {1, 0, 0});
    test({mesh, Ms_field, m({0, 1, 0})}, {0, 1, 0});
    test({mesh, Ms_field, m({0, 0, 1})}, {0, 0, 1});

    // Ms.field with zeros
    auto Ms_field_with_zeros = af::constant(0, dims_scalar(mesh), f64);
    Ms_field_with_zeros(0, 0, 0) = Ms;
    test({mesh, Ms_field_with_zeros, m({1, 0, 0})}, {1, 0, 0});
    test({mesh, Ms_field_with_zeros, m({0, 1, 0})}, {0, 1, 0});
    test({mesh, Ms_field_with_zeros, m({0, 0, 1})}, {0, 0, 1});

    // This should work:
    // auto res = State(mesh, Ms_field_with_zeros, m({0, 0, 1})).mean_m();
    // EXPECT_THAT(res, {0, 0, 1});
    // EXPECT_THAT(res, std::array<double, 3>({0, 0, 1}));
}
int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
