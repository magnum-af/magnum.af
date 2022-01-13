#include "state.hpp"
#include "util/util.hpp"
#include <cmath> // std::isnan
#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <iostream>
#include <vector>

using namespace magnumaf;
TEST(State, Init_Ms_field) {
    int nx = 1, ny = 1, nz = 1;
    auto m = af::constant(0, nx, ny, nz, 3, f64);
    m(af::span, af::span, af::span, 0) = 1;
    auto Ms_scalar = af::constant(1, nx, ny, nz, 1, f64);
    State state_scalar(Mesh(1, 1, 1, 1, 1, 1), Ms_scalar, m);
    EXPECT_EQ(state_scalar.Ms_field.dims(3), 1);

    auto Ms_vector = af::constant(1, nx, ny, nz, 1, f64);
    State state_legacy(Mesh(1, 1, 1, 1, 1, 1), Ms_vector, m);
    EXPECT_EQ(state_legacy.Ms_field.dims(3), 1);

    // To test wrapping we use python, elso somting like this:
    // auto Ms_wrapping = af::constant(1, nx, ny, nz, 1, f64);
    // State state_wrapping(Mesh(1,1,1,1,1,1), (long int) Ms_scalar.get(), (long
    // int) m.get()); EXPECT_EQ(state_wrapping.Ms_field.dims(3), 3);
}

TEST(State, Ms_is_quiet_Nan_if_Msfield_is_set) {
    int nx = 1, ny = 1, nz = 1;
    auto m = af::constant(0, nx, ny, nz, 3, f64);
    auto Ms_scalar = af::constant(1, nx, ny, nz, 1, f64);
    State state(Mesh(1, 1, 1, 1, 1, 1), Ms_scalar, m);
    EXPECT_EQ(std::isnan(state.Ms), true);
    EXPECT_EQ(std::isnan(2 * state.Ms), true);
}

TEST(State, normalize_m) {
    Mesh mesh(1, 1, 1, 1e-9, 1e-9, 1e-9);
    double Ms = 1;
    // af::array m = af::constant(std::sqrt(3), 1, 1, 1, 3, f64); // TODO check why this prints custom warning,
    // numerical?
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

    auto test = [](const State& state, std::array<double, 3> res) { EXPECT_THAT(state.mean_m(), res); };

    // vec_as_array method:
    auto m = [dims = mesh::dims_s(mesh)](std::array<double, 3> vec) {
        return af::tile(af::array(af::dim4(1, 1, 1, 3), vec.data()), dims);
    };

    // scalar Ms
    {
        test(State{mesh, Ms, m({1, 0, 0})}, {1, 0, 0});
        test(State{mesh, Ms, m({0, 1, 0})}, {0, 1, 0});
        test(State{mesh, Ms, m({0, 0, 1})}, {0, 0, 1});
    }

    // Ms.field
    {
        auto Ms_field = af::constant(Ms, mesh::dims_s(mesh), f64);
        test(State{mesh, Ms_field, m({1, 0, 0})}, {1, 0, 0});
        test(State{mesh, Ms_field, m({0, 1, 0})}, {0, 1, 0});
        test(State{mesh, Ms_field, m({0, 0, 1})}, {0, 0, 1});
    }

    // Ms.field with zeros
    {
        auto Ms_field_with_zeros = af::constant(0, mesh::dims_s(mesh), f64);
        Ms_field_with_zeros(0, 0, 0) = Ms;
        test(State{mesh, Ms_field_with_zeros, m({1, 0, 0})}, {1, 0, 0});
        test(State{mesh, Ms_field_with_zeros, m({0, 1, 0})}, {0, 1, 0});
        test(State{mesh, Ms_field_with_zeros, m({0, 0, 1})}, {0, 0, 1});
    }

    // Ms.field and m with zeros
    {
        auto Ms_field_with_zeros = af::constant(0, mesh::dims_s(mesh), f64);
        Ms_field_with_zeros(0, 0, 0) = Ms;

        {
            auto m0 = af::constant(0, mesh::dims_v(mesh), f64);
            m0(0, af::span, af::span, 0) = 1;
            test(State{mesh, Ms_field_with_zeros, m0}, {1, 0, 0});
        }

        {
            auto m0 = af::constant(0, mesh::dims_v(mesh), f64);
            m0(0, af::span, af::span, 1) = 1;
            test(State{mesh, Ms_field_with_zeros, m0}, {0, 1, 0});
        }

        {
            auto m0 = af::constant(0, mesh::dims_v(mesh), f64);
            m0(0, af::span, af::span, 2) = 1;
            test(State{mesh, Ms_field_with_zeros, m0}, {0, 0, 1});
        }
    }

    // m with zeros
    {
        {
            auto m0 = af::constant(0, mesh::dims_v(mesh), f64);
            m0(0, af::span, af::span, 0) = 1;
            test(State{mesh, Ms, m0}, {1, 0, 0});
        }

        {
            auto m0 = af::constant(0, mesh::dims_v(mesh), f64);
            m0(0, af::span, af::span, 1) = 1;
            test(State{mesh, Ms, m0}, {0, 1, 0});
        }

        {
            auto m0 = af::constant(0, mesh::dims_v(mesh), f64);
            m0(0, af::span, af::span, 2) = 1;
            test(State{mesh, Ms, m0}, {0, 0, 1});
        }
    }

    // This should work:
    // auto res = State(mesh, Ms_field_with_zeros, m({0, 0, 1})).mean_m();
    // EXPECT_THAT(res, {0, 0, 1});
    // EXPECT_THAT(res, std::array<double, 3>({0, 0, 1}));
}

TEST(State, mean_M) {
    const auto get_mean_M = [](auto const& Ms_field, af::array const& m0) {
        auto nx = m0.dims(0);
        auto ny = m0.dims(0);
        auto nz = m0.dims(0);
        Mesh mesh(nx, ny, nz, 1e-9, 1e-9, 1e-9);
        auto state = State(mesh, Ms_field, m0);
        auto mean_M = state.mean_M();
        return mean_M;
    };
    double Ms = 1e6;
    const af::dim4 dim_s(4, 1, 1, 1);
    const af::dim4 dim_v(4, 1, 1, 3);
    {
        auto m0 = af::constant(0.0, dim_v, f64);
        m0(af::span, af::span, af::span, 0) = 1;
        auto mean_M = get_mean_M(af::constant(Ms, dim_s, f64), m0);
        EXPECT_THAT(mean_M, testing::ElementsAre(Ms, 0, 0));
    }
    // Testing zero-values m0
    {
        auto m0 = af::constant(0.0, dim_v, f64);
        m0(0, 0, 0, 0) = 1;
        m0(1, 0, 0, 0) = 1;
        auto mean_M = get_mean_M(af::constant(Ms, dim_s, f64), m0);
        EXPECT_THAT(mean_M, testing::ElementsAre(Ms, 0, 0));
    }
    // Testing zero-values Ms_field
    {
        auto Ms_field = af::constant(0.0, dim_s, f64);
        Ms_field(af::seq(0, 1), 0, 0) = Ms;
        auto m0 = af::constant(0.0, dim_v, f64);
        m0(af::span, af::span, af::span, 0) = 1;
        auto mean_M = get_mean_M(Ms_field, m0);
        EXPECT_THAT(mean_M, testing::ElementsAre(Ms / 2, 0, 0));
    }
    // Testing <Ms>=0
    {
        auto m0 = af::constant(0.0, dim_v, f64);
        m0(af::seq(0, 1), af::span, af::span, 1) = 1;
        m0(af::seq(2, 3), af::span, af::span, 1) = -1;
        auto mean_M = get_mean_M(af::constant(Ms, dim_s, f64), m0);
        EXPECT_THAT(mean_M, testing::ElementsAre(0, 0, 0));
    }
    // Testing scalar Ms (i.e. no Ms_field)
    {
        auto m0 = af::constant(0.0, dim_v, f64);
        m0(af::span, af::span, af::span, 1) = 1;
        auto mean_M = get_mean_M(Ms, m0);
        EXPECT_THAT(mean_M, testing::ElementsAre(0, Ms, 0));
    }
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
