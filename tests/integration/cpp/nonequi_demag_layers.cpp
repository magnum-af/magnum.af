#include <gtest/gtest.h>
#include "../../../src/mesh.cpp"
#include "../../../src/constants.hpp"
#include "../../../src/material.cpp"
#include "../../../src/state.cpp"
#include "../../../src/vtk_IO.cpp"
#include "../../../src/func.cpp"
#include "../../../src/misc.cpp"
#include "../../../src/llg_terms/micro_demag_nonequi.cpp"
#include "../../../src/llg_terms/micro_demag.cpp"
 
TEST(NonEquiDemagSP4LayerTest, n) {
    // Compare SP4 layer (with ↑ → → → → ↑ magnetization) once with nz = 3 equidistant and nz=2 non-equidistant discretization
    const double x=5.e-7, y=1.25e-7, z=3.e-9;
    const int nx = 100, ny=25, nz=3, nz_ne=2;

    // equi
    Mesh mesh_ed(nx,ny,nz,x/nx,y/ny,z/nz);
    Material material_ed = Material();
    material_ed.ms    = 8e5;

    af::array m = af::constant(0.0,mesh_ed.n0,mesh_ed.n1,mesh_ed.n2,3,f64);
    m(af::seq(1,af::end-1),af::span,af::span,0) = af::constant(1.0,mesh_ed.n0-2,mesh_ed.n1,mesh_ed.n2,1,f64);
    m(0,af::span,af::span,1 ) = af::constant(1.0,1,mesh_ed.n1,mesh_ed.n2,1,f64);
    m(-1,af::span,af::span,1) = af::constant(1.0,1,mesh_ed.n1,mesh_ed.n2,1,f64);
    //causes error of 1e-2//m(af::span, af::span, 0, af::span) = 0;

    State state_ed(mesh_ed, material_ed, m);
    DemagField demag_ed = DemagField(mesh_ed, material_ed, false, false, 1);

    // nonequi
    std::vector<double> z_spacing = {z/nz, 2 * z/nz};
    Mesh mesh_ne(nx, ny, nz_ne, x/nx, y/ny, z/nz_ne);
    Material material_ne = Material();
    material_ne.ms    = 8e5;
    af::array m2 = af::constant(0.0,mesh_ne.n0,mesh_ne.n1,mesh_ne.n2,3,f64);
    m2(af::seq(1,af::end-1),af::span,af::span,0) = af::constant(1.0,mesh_ne.n0-2,mesh_ne.n1,mesh_ne.n2,1,f64);
    m2(0,af::span,af::span,1 ) = af::constant(1.0,1,mesh_ne.n1,mesh_ne.n2,1,f64);
    m2(-1,af::span,af::span,1) = af::constant(1.0,1,mesh_ne.n1,mesh_ne.n2,1,f64);
    //causes error of 1e-2//m2(af::span, af::span, 0, af::span) = 0;

    State state_ne(mesh_ne, material_ne, m2);
    NonEquiDemagField demag_ne = NonEquiDemagField(mesh_ne, z_spacing, false, false, 1);

    EXPECT_NEAR( max_abs_diff(constants::mu0 * demag_ed.h(state_ed)(af::span, af::span, 0, af::span), constants::mu0 * demag_ne.h(state_ne)(af::span, af::span, 0, af::span)), 0, 1e-08);
    EXPECT_NEAR(mean_abs_diff(constants::mu0 * demag_ed.h(state_ed)(af::span, af::span, 0, af::span), constants::mu0 * demag_ne.h(state_ne)(af::span, af::span, 0, af::span)), 0, 1e-08);

    //vti_writer_micro(state_ed.m, mesh_ed ,"m_equi");
    //vti_writer_micro(state_ne.m, mesh_ne ,"m_nonequi");
    //vti_writer_micro(demag_ed.h(state_ed), mesh_ed ,"h_equi");
    //vti_writer_micro(demag_ne.h(state_ne), mesh_ne ,"h_nonequi");
}
TEST(NonEquiDemagHomogenLayerTest, n) {
    // Compare SP4 layer with homogenuous z-magnetization once with nz = 3 equidistant and nz=2 non-equidistant discretization
    const double x=5.e-7, y=1.25e-7, z=3.e-9;
    const int nx = 100, ny=25, nz=3, nz_ne=2;

    // equi
    Mesh mesh_ed(nx,ny,nz,x/nx,y/ny,z/nz);
    Material material_ed = Material();
    material_ed.ms    = 8e5;

    af::array m = af::constant(0.0,mesh_ed.n0,mesh_ed.n1,mesh_ed.n2,3,f64);
    m(af::span, af::span, af::span, 2) = 1;

    State state_ed(mesh_ed, material_ed, m);
    DemagField demag_ed = DemagField(mesh_ed, material_ed, false, false, 1);

    // nonequi
    std::vector<double> z_spacing = {z/nz, 2 * z/nz};
    Mesh mesh_ne(nx, ny, nz_ne, x/nx, y/ny, z/nz_ne);
    Material material_ne = Material();
    material_ne.ms    = 8e5;
    af::array m2 = af::constant(0.0,mesh_ne.n0,mesh_ne.n1,mesh_ne.n2,3,f64);
    m2(af::span, af::span, af::span, 2) = 1;

    State state_ne(mesh_ne, material_ne, m2);
    NonEquiDemagField demag_ne = NonEquiDemagField(mesh_ne, z_spacing, false, false, 1);

    EXPECT_NEAR(max_abs_diff(demag_ed.h(state_ed)(af::span, af::span, 0, af::span), demag_ne.h(state_ne)(af::span, af::span, 0, af::span)) * constants::mu0, 0, 1e-8);
    EXPECT_NEAR(mean_abs_diff(constants::mu0 * demag_ed.h(state_ed)(af::span, af::span, 0, af::span), constants::mu0 * demag_ne.h(state_ne)(af::span, af::span, 0, af::span)), 0, 1e-08);

    //vti_writer_micro(state_ed.m, mesh_ed ,"m_equi");
    //vti_writer_micro(state_ne.m, mesh_ne ,"m_nonequi");
    //vti_writer_micro(demag_ed.h(state_ed), mesh_ed ,"h_equi");
    //vti_writer_micro(demag_ne.h(state_ne), mesh_ne ,"h_nonequi");
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
