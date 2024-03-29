#include "field_terms/micro/demag_field.hpp"
#include "field_terms/nonequi/nonequi_demag_field.hpp"
#include "util/util.hpp"
#include <gtest/gtest.h>

using namespace magnumaf;

TEST(NonequiDemagField, energy_homogenuous_cube) {
    // Compare SP4 layer with homogeneous z-magnetization once with nz = 3
    // equidistant and nz=2 non-equidistant discretization
    const double dx = 1e-9, dy = 1e-9, dz = 1e-9;
    const int nx = 4, ny = 4, nz = 4;
    const double Ms = 8e5;

    std::vector<double> z_spacing;
    for (int i = 0; i < nz; i++) {
        if (i % 2 == 0) {
            z_spacing.push_back(0.5 * dz);
        } else {
            z_spacing.push_back(1.5 * dz);
        }
    }
    NonequiMesh mesh_ne(nx, ny, dx, dy, z_spacing);
    af::array m2 = af::constant(0.0, nemesh::dims_v(mesh_ne), f64);
    m2(af::span, af::span, af::span, 2) = 1;

    State state(m2, Ms, false, true);
    NonequiDemagField demag = NonequiDemagField(mesh_ne, false, false, 1);

    double E_analytic = 1. / 6. * nx * pow(dx, 3) * pow(state.Ms, 2) * constants::mu0;

    af::array h_demag = demag.H_in_Apm(state);
    // testing state.Ms
    EXPECT_NEAR(demag.Energy_in_J(state), E_analytic, 3.7e-19);
    EXPECT_NEAR(demag.Energy_in_J(state, h_demag), E_analytic, 3.7e-19);

    // testing state.Ms_field switch
    state.Ms_field = af::constant(8e5, nemesh::dims_s(mesh_ne), f64);
    h_demag = demag.H_in_Apm(state);
    EXPECT_NEAR(demag.Energy_in_J(state), E_analytic, 3.7e-20); // 10x10x10 := 3.7e-19
    EXPECT_NEAR(demag.Energy_in_J(state, h_demag), E_analytic, 3.7e-20);
}

TEST(NonequiDemagField, EnergyTest) {
    // Compare SP4 layer with random z-magnetization once with nz = 3
    // equidistant and nz=2 non-equidistant discretization NOTE: this test is
    // only sensitive to i_source >= i_target (in NonequiDemagField::H_in_Apm() method)
    const double x = 5.e-7, y = 1.25e-7, z = 3.e-9;
    const int nx = 10, ny = 10, nz = 3;

    af::randomEngine rand_engine = util::rand_engine_current_time();
    af::array random_1 = af::randu(af::dim4(nx, ny, 1, 3), f64, rand_engine);
    af::array random_2 = af::randu(af::dim4(nx, ny, 1, 3), f64, rand_engine);

    // equi
    Mesh mesh_ed(nx, ny, nz, x / nx, y / ny, z / nz);

    af::array m = af::constant(0.0, mesh_ed.nx, mesh_ed.ny, mesh_ed.nz, 3, f64);
    m(af::span, af::span, 0, af::span) = random_1;
    m(af::span, af::span, 1, af::span) = random_2;
    m(af::span, af::span, 2, af::span) = random_2;

    State state_ed(mesh_ed, 8e5, m, false, true);
    DemagField demag_ed = DemagField(mesh_ed, false, false, 1);

    // nonequi
    std::vector<double> z_spacing = {z / nz, 2 * z / nz};
    NonequiMesh mesh_ne(nx, ny, x / nx, y / ny, z_spacing);
    af::array m2 = af::constant(0.0, nemesh::dims_v(mesh_ne), f64);
    m2(af::span, af::span, 0, af::span) = random_1;
    m2(af::span, af::span, 1, af::span) = random_2;

    // TODO this is not throwing an error but segfaults. No compiler warnign
    // even though constructor should not exist: //State state_ne(mesh_ne,
    // (af::array) m2, false, true);
    State state_ne(m2, af::constant(8e5, nemesh::dims_s(mesh_ne), f64), false, true);
    NonequiDemagField demag_ne = NonequiDemagField(mesh_ne, false, false, 1);

    // testing E(state)
    EXPECT_NEAR(demag_ed.Energy_in_J(state_ed), demag_ne.Energy_in_J(state_ne),
                2e-24); // Note: this is for opencl, cpu and cuda achieve 3e-27

    // testing E(state, heff)
    af::array h_ed = demag_ed.H_in_Apm(state_ed);
    af::array h_ne = demag_ne.H_in_Apm(state_ne);
    EXPECT_NEAR(demag_ed.Energy_in_J(state_ed, h_ed), demag_ne.Energy_in_J(state_ne, h_ne),
                2e-24); // Note: this is for opencl, cpu and cuda achieve 3e-27
}

TEST(NonequiDemagField, RandomMagnetizationHeffTest) {
    // Compare SP4 layer with random z-magnetization once with nz = 3
    // equidistant and nz=2 non-equidistant discretization NOTE: this test is
    // only sensitive to i_source >= i_target (in NonequiDemagField::H_in_Apm() method)
    const double x = 5.e-7, y = 1.25e-7, z = 3.e-9;
    const int nx = 10, ny = 10, nz = 3;

    af::randomEngine rand_engine = util::rand_engine_current_time();
    af::array random_1 = af::randu(af::dim4(nx, ny, 1, 3), f64, rand_engine);
    af::array random_2 = af::randu(af::dim4(nx, ny, 1, 3), f64, rand_engine);

    // equi
    Mesh mesh_ed(nx, ny, nz, x / nx, y / ny, z / nz);

    af::array m = af::constant(0.0, mesh_ed.nx, mesh_ed.ny, mesh_ed.nz, 3, f64);
    m(af::span, af::span, 0, af::span) = random_1;
    m(af::span, af::span, 1, af::span) = random_2;
    m(af::span, af::span, 2, af::span) = random_2;

    State state_ed(mesh_ed, 8e5, m, false, true);
    DemagField demag_ed = DemagField(mesh_ed, false, false, 1);

    // nonequi
    std::vector<double> z_spacing = {z / nz, 2 * z / nz};
    NonequiMesh mesh_ne(nx, ny, x / nx, y / ny, z_spacing);
    af::array m2 = af::constant(0.0, nemesh::dims_v(mesh_ne), f64);
    m2(af::span, af::span, 0, af::span) = random_1;
    m2(af::span, af::span, 1, af::span) = random_2;

    State state_ne(m2, af::constant(8e5, nemesh::dims_s(mesh_ne), f64), false, true);
    NonequiDemagField demag_ne = NonequiDemagField(mesh_ne, false, false, 1);

    af::array demag_ed_h = demag_ed.H_in_Apm(state_ed)(af::span, af::span, 0, af::span);
    af::array demag_ne_h = demag_ne.H_in_Apm(state_ne)(af::span, af::span, 0, af::span);
    // af::print("demag_ed_h", demag_ed_h);
    // af::print("demag_ne_h", demag_ne_h);
    EXPECT_NEAR(util::max_abs_diff(demag_ed_h, demag_ne_h), 0,
                0.04); // Note: this is for opencl, cpu and cuda achieve 0.007
    EXPECT_NEAR(util::mean_abs_diff(demag_ed_h, demag_ne_h), 0,
                0.009); // Note: this is for opencl, cpu and cuda achieve
}

TEST(NonequiDemagField, RandomMagnetizationSwappedZindexHeffTest) {
    // Same as above but testing else{} in nedemag.H_in_Apm(state) method
    // NOTE: this test is only sensitive to i_source < i_target (in
    // NonequiDemagField::H_in_Apm() method)
    const double x = 5.e-7, y = 1.25e-7, z = 3.e-9;
    const int nx = 10, ny = 10, nz = 3;

    af::randomEngine rand_engine = util::rand_engine_current_time();
    af::array random_1 = af::randu(af::dim4(nx, ny, 1, 3), f64, rand_engine);
    af::array random_2 = af::randu(af::dim4(nx, ny, 1, 3), f64, rand_engine);

    // equi
    Mesh mesh_ed(nx, ny, nz, x / nx, y / ny, z / nz);

    af::array m = af::constant(0.0, mesh_ed.nx, mesh_ed.ny, mesh_ed.nz, 3, f64);
    m(af::span, af::span, 0, af::span) = random_2;
    m(af::span, af::span, 1, af::span) = random_2;
    m(af::span, af::span, 2, af::span) = random_1;

    State state_ed(mesh_ed, 8e5, m, false, true);
    DemagField demag_ed = DemagField(mesh_ed, false, false, 1);

    // nonequi
    std::vector<double> z_spacing = {2 * z / nz, z / nz};
    NonequiMesh mesh_ne(nx, ny, x / nx, y / ny, z_spacing);
    af::array m2 = af::constant(0.0, nemesh::dims_v(mesh_ne), f64);
    m2(af::span, af::span, 0, af::span) = random_2;
    m2(af::span, af::span, 1, af::span) = random_1;

    State state_ne(m2, af::constant(8e5, nemesh::dims_s(mesh_ne), f64), false, true);
    NonequiDemagField demag_ne = NonequiDemagField(mesh_ne, false, false, 1);

    af::array demag_ed_h = demag_ed.H_in_Apm(state_ed)(af::span, af::span, 2, af::span);
    af::array demag_ne_h = demag_ne.H_in_Apm(state_ne)(af::span, af::span, 1, af::span);
    // af::print("demag_ed_h", demag_ed_h);
    // af::print("demag_ne_h", demag_ne_h);
    EXPECT_NEAR(util::max_abs_diff(demag_ed_h, demag_ne_h), 0,
                0.04); // Note: this is for opencl, cpu and cuda achieve 0.007
    EXPECT_NEAR(util::mean_abs_diff(demag_ed_h, demag_ne_h), 0,
                0.009); // Note: this is for opencl, cpu and cuda achieve 0.0009
}

TEST(NonequiDemagField, RandomMagnetizationWithZeroLayerHeffTest) {
    // Compare SP4 layer with homogeneous z-magnetization once with nz = 3
    // equidistant and nz=2 non-equidistant discretization
    const double x = 5.e-7, y = 1.25e-7, z = 3.e-9;
    const int nx = 10, ny = 8, nz = 3;

    af::randomEngine rand_engine = util::rand_engine_current_time();
    af::array random = af::randu(af::dim4(nx, ny, 1, 3), f64, rand_engine);

    // equi
    Mesh mesh_ed(nx, ny, nz, x / nx, y / ny, z / nz);

    af::array m = af::constant(0.0, mesh_ed.nx, mesh_ed.ny, mesh_ed.nz, 3, f64);
    m(af::span, af::span, 1, af::span) = random;
    m(af::span, af::span, 2, af::span) = random;

    State state_ed(mesh_ed, 8e5, m, false,
                   true); // Note: this implicitly sets Ms_filed internally
    DemagField demag_ed = DemagField(mesh_ed, false, false, 1);

    // nonequi
    std::vector<double> z_spacing = {z / nz, 2 * z / nz};
    NonequiMesh mesh_ne(nx, ny, x / nx, y / ny, z_spacing);
    af::array m2 = af::constant(0.0, nemesh::dims_v(mesh_ne), f64);
    m2(af::span, af::span, 1, af::span) = random;

    // af::array Ms_field = af::constant(8e5, nemesh::dims_s(mesh_ne), f64);
    // std::cout << "test" << std::endl;
    // State state_ne(mesh_ne, Ms_field, m2, false, true);
    // std::cout << "test" << std::endl;
    State state_ne(m2, 8e5, false, true);
    NonequiDemagField demag_ne = NonequiDemagField(mesh_ne, false, false, 1);

    af::array demag_ed_h = demag_ed.H_in_Apm(state_ed)(af::span, af::span, 0, af::span);
    af::array demag_ne_h = demag_ne.H_in_Apm(state_ne)(af::span, af::span, 0, af::span);
    EXPECT_NEAR(util::max_abs_diff(demag_ed_h, demag_ne_h), 0,
                0.004); // cpu and cuda 0.001; 100x25x3: 0.01
    EXPECT_NEAR(util::mean_abs_diff(demag_ed_h, demag_ne_h), 0,
                0.0009); // cpu and cuda 0.0003; 100x25x3:  0.003
}

TEST(NonequiDemagField, UMagnetizationHeffTest) {
    // Compare SP4 layer (with ↑ → → → → ↑ magnetization) once with nz = 3
    // equidistant and nz=2 non-equidistant discretization
    const double x = 5.e-7, y = 1.25e-7, z = 3.e-9;
    const int nx = 10, ny = 8, nz = 3;

    // equi
    Mesh mesh_ed(nx, ny, nz, x / nx, y / ny, z / nz);

    af::array m = af::constant(0.0, mesh_ed.nx, mesh_ed.ny, mesh_ed.nz, 3, f64);
    m(af::seq(1, af::end - 1), af::span, af::span, 0) =
        af::constant(1.0, mesh_ed.nx - 2, mesh_ed.ny, mesh_ed.nz, 1, f64);
    m(0, af::span, af::span, 1) = af::constant(1.0, 1, mesh_ed.ny, mesh_ed.nz, 1, f64);
    m(-1, af::span, af::span, 1) = af::constant(1.0, 1, mesh_ed.ny, mesh_ed.nz, 1, f64);
    m(af::span, af::span, 0, af::span) = 0;

    State state_ed(mesh_ed, 8e5, m, false, true);
    DemagField demag_ed = DemagField(mesh_ed, false, false, 1);

    // nonequi
    std::vector<double> z_spacing = {z / nz, 2 * z / nz};
    NonequiMesh mesh_ne(nx, ny, x / nx, y / ny, z_spacing);
    af::array m2 = af::constant(0.0, nemesh::dims_v(mesh_ne), f64);
    m2(af::seq(1, af::end - 1), af::span, af::span, 0) =
        af::constant(1.0, mesh_ne.nx - 2, mesh_ne.ny, mesh_ne.nz, 1, f64);
    m2(0, af::span, af::span, 1) = af::constant(1.0, 1, mesh_ne.ny, mesh_ne.nz, 1, f64);
    m2(-1, af::span, af::span, 1) = af::constant(1.0, 1, mesh_ne.ny, mesh_ne.nz, 1, f64);
    m2(af::span, af::span, 0, af::span) = 0;

    State state_ne(m2, 8e5, false, true);
    state_ne.Ms_field = af::constant(8e5, nemesh::dims_s(mesh_ne), f64);
    NonequiDemagField demag_ne = NonequiDemagField(mesh_ne, false, false, 1);

    af::array demag_ed_h = demag_ed.H_in_Apm(state_ed)(af::span, af::span, 0, af::span);
    af::array demag_ne_h = demag_ne.H_in_Apm(state_ne)(af::span, af::span, 0, af::span);
    EXPECT_NEAR(util::max_abs_diff(demag_ed_h, demag_ne_h), 0,
                0.002); // 0.0001; 100x25x3: 0.01, Note: this is for opencl, cpu and
                        // cuda achieve 0.007
    EXPECT_NEAR(util::mean_abs_diff(demag_ed_h, demag_ne_h), 0,
                0.0002); // 0.000008; 100x25x3: 0.0008, Note: this is for opencl, cpu
                         // and cuda achieve 0.0006
}

TEST(NonequiDemagField, HomogenuousMagnetizationHeffTest) {
    // Compare SP4 layer with homogeneous z-magnetization once with nz = 3
    // equidistant and nz=2 non-equidistant discretization
    const double x = 5.e-7, y = 1.25e-7, z = 3.e-9;
    const int nx = 10, ny = 8, nz = 3;

    // equi
    Mesh mesh_ed(nx, ny, nz, x / nx, y / ny, z / nz);

    af::array m = af::constant(0.0, mesh_ed.nx, mesh_ed.ny, mesh_ed.nz, 3, f64);
    m(af::span, af::span, af::span, 2) = 1;

    State state_ed(mesh_ed, 8e5, m, false, true);
    DemagField demag_ed = DemagField(mesh_ed, false, false, 1);

    // nonequi
    std::vector<double> z_spacing = {z / nz, 2 * z / nz};
    NonequiMesh mesh_ne(nx, ny, x / nx, y / ny, z_spacing);
    af::array m2 = af::constant(0.0, nemesh::dims_v(mesh_ne), f64);
    m2(af::span, af::span, af::span, 2) = 1;

    State state_ne(m2, 8e5, false, true);
    state_ne.Ms_field = af::constant(8e5, nemesh::dims_s(mesh_ne), f64);
    NonequiDemagField demag_ne = NonequiDemagField(mesh_ne, false, false, 1);

    af::array demag_ed_h = demag_ed.H_in_Apm(state_ed)(af::span, af::span, 0, af::span);
    af::array demag_ne_h = demag_ne.H_in_Apm(state_ne)(af::span, af::span, 0, af::span);
    EXPECT_NEAR(util::max_abs_diff(demag_ed_h, demag_ne_h), 0,
                0.04); // 0.0004; 100x25x3: 0.04 Note: this is for opencl, cpu and
                       // cuda achieve 0.007
    EXPECT_NEAR(util::mean_abs_diff(demag_ed_h, demag_ne_h), 0,
                0.01); // 0.0001; 100x25x3: 0.01 Note: this is for opencl, cpu and
                       // cuda achieve 0.0003
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
