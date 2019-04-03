#include "magnum_af.hpp"
#include "llg_terms/micro_demag.hpp"
#include "llg_terms/micro_nonequi_demag.hpp"

int main(int argc, char** argv)
{
    const int run_test = 3; // specify test no., 0 is all

    //Test 1: Testing heff for 
    if (run_test == 0 || run_test == 1)
    {
        const double x=5.e-7, y=1.25e-7, z=3.e-9;
        //const int nx = 100, ny=25 ,nz=1;
        const int nx = 3, ny=1 ,nz=2;
        const int nz_nonequi = 1;
        
        //Mesh mesh(nx,ny,nz,x/nx,y/ny,z/nz);
        Material material = Material();
        material.ms    = 8e5;
        material.A     = 1.3e-11;
        material.alpha = 1;

        Mesh mesh_full(nx, ny, nz, x/nx,y/ny,z/nz);
        Mesh mesh_nonequi(nx, ny, nz_nonequi, x/nx,y/ny,z/nz);


        af::array m_nonequi = af::constant(0, nx, ny, nz_nonequi, 3, f64);

        m_nonequi(af::seq(1,af::end-1),af::span,af::span,0) = af::constant(1.0,nx-2,ny,nz_nonequi,1,f64);
        m_nonequi(0,af::span,af::span,1 ) = af::constant(1.0,1,ny,nz_nonequi,1,f64);
        m_nonequi(-1,af::span,af::span,1) = af::constant(1.0,1,ny,nz_nonequi,1,f64);

        State state_nonequi(mesh_nonequi, material, m_nonequi);

        af::array m_full = af::constant(0, nx, ny, nz, 3, f64);
        m_full(af::seq(1,af::end-1),af::span, 0, 0) = af::constant(1.0,nx-2,ny, 1,1,f64);
        m_full( 0, af::span, 0, 1) = af::constant(1.0, 1, ny, 1, 1, f64);
        m_full(-1, af::span, 0, 1) = af::constant(1.0, 1, ny, 1, 1, f64);
        State state_full(mesh_full,material, m_full);

        DemagField demag = DemagField(mesh_full, material, true, false, 1);

        NonEquiDemagField nonequi_demag = NonEquiDemagField(mesh_nonequi, material, true, false, 1);

        //af::print("full h", demag.h(state_full));
        af::print("full h", demag.h(state_full)(af::span, af::span, 1, af::span));

        af::print("nonequi_demag H", nonequi_demag.h(state_nonequi));

        //af::print("demag.Nfft", demag.Nfft);
        //af::print("nonequi_demag.Nfft", nonequi_demag.Nfft);
        //unsigned int zero_if_equal_var = zero_if_equal(demag.Nfft, nonequi_demag.Nfft);
        
        abs_diff_lt_precision(demag.h(state_full)(af::span, af::span, 1, af::span), nonequi_demag.h(state_nonequi));
        //if (abs_diff_lt_precision(demag.h(state_full)(af::span, af::span, 1, af::span), 2*nonequi_demag.h(state_nonequi))) std::cout << "test true" << std::endl;
        //else std::cout << "test false" << std::endl;
        //zero_if_equal(demag.Nfft, nonequi_demag.Nfft);
        //zero_if_equal(demag.Nfft(af::span, af::span, 2, af::span), nonequi_demag.Nfft);
    }

    //Test 2: Testing heff with extendend newell 
    //NOTE: requires j2 = 1 in /micro_nonequi_demag.cpp:272
    if (run_test == 0 || run_test == 2)
    {
        const double x=5.e-7, y=1.25e-7, z=3.e-9;
        //const int nx = 100, ny=25 ,nz=1;
        const int nx = 3, ny=1 ,nz=2;
        const int nz_nonequi = 1;
        
        //Mesh mesh(nx,ny,nz,x/nx,y/ny,z/nz);
        Material material = Material();
        material.ms    = 8e5;
        material.A     = 1.3e-11;
        material.alpha = 1;

        Mesh mesh_full(nx, ny, nz, x/nx,y/ny,z/nz);
        Mesh mesh_nonequi(nx, ny, nz_nonequi, x/nx,y/ny,z/nz);


        af::array m_nonequi = af::constant(0, nx, ny, nz_nonequi, 3, f64);

        m_nonequi(af::seq(1,af::end-1),af::span,af::span,0) = af::constant(1.0,nx-2,ny,nz_nonequi,1,f64);
        m_nonequi(0,af::span,af::span,1 ) = af::constant(1.0,1,ny,nz_nonequi,1,f64);
        m_nonequi(-1,af::span,af::span,1) = af::constant(1.0,1,ny,nz_nonequi,1,f64);

        State state_nonequi(mesh_nonequi, material, m_nonequi);

        af::array m_full = af::constant(0, nx, ny, nz, 3, f64);
        m_full(af::seq(1,af::end-1),af::span, 0, 0) = af::constant(1.0,nx-2,ny, 1,1,f64);
        m_full( 0, af::span, 0, 1) = af::constant(1.0, 1, ny, 1, 1, f64);
        m_full(-1, af::span, 0, 1) = af::constant(1.0, 1, ny, 1, 1, f64);
        State state_full(mesh_full,material, m_full);

        DemagField demag = DemagField(mesh_full, material, true, false, 1);

        NonEquiDemagField nonequi_demag = NonEquiDemagField(mesh_nonequi, material, true, false, 1);

        af::print("full h", demag.h(state_full)(af::span, af::span, 1, af::span));

        af::print("nonequi_demag H", nonequi_demag.h(state_nonequi));

        abs_diff_lt_precision(demag.h(state_full)(af::span, af::span, 1, af::span), nonequi_demag.h(state_nonequi), 3.9e-8);
        abs_diff_lt_precision(demag.h(state_full)(af::span, af::span, 1, af::span), nonequi_demag.h(state_nonequi), 4e-8);
        rel_diff_lt_precision(demag.h(state_full)(af::span, af::span, 1, af::span), nonequi_demag.h(state_nonequi), 1e-03);
        rel_diff_lt_precision(demag.h(state_full)(af::span, af::span, 1, af::span), nonequi_demag.h(state_nonequi), 2e-03);

        //test
        //rel_diff_lt_precision(af::constant(1., 1, f64), af::constant(1.1, 1, f64), 0.1);
        //rel_diff_lt_precision(af::constant(1., 1, f64), af::constant(1.1, 1, f64), 0.099);
        //rel_diff_lt_precision(af::constant(1., 1, f64), af::constant(1.1, 1, f64), 0.09);
        //rel_diff_lt_precision(af::constant(-1., 1, f64), af::constant(-1.1, 1, f64), 0.1);
        //rel_diff_lt_precision(af::constant(1., 1, f64), af::constant(1.1, 1, f64), 0.1);
        //rel_diff_lt_precision(af::constant(10., 1, f64), af::constant(10.1, 1, f64), 0.1);
    }
    //Test 3: Testing heff with extendend newell  for sp4, not yet working
    //Note: Negative check with using old newell sucessfull
    if (run_test == 0 || run_test == 3)
    {
        const double x=5.e-7, y=1.25e-7, z=3.e-9;
        const int nx = 100, ny=25 ,nz=1;
        //const int nx = 3, ny=1 ,nz=2;
        //const int nx = 3, ny=1 ,nz=1;
        //const int nx = 3, ny=3 ,nz=3;
        
        Mesh mesh(nx,ny,nz,x/nx,y/ny,z/nz);
        Material material = Material();
        material.ms    = 8e5;
        material.A     = 1.3e-11;
        material.alpha = 1;

        State state_nonequi(mesh, material, mesh.init_sp4());

        State state_full(mesh, material, mesh.init_sp4());

        DemagField demag = DemagField(mesh, material, true, false, 1);

        NonEquiDemagField nonequi_demag = NonEquiDemagField(mesh, material, true, false, 1);

        //af::print("full h", demag.h(state_full));

        //af::print("nonequi_demag H", nonequi_demag.h(state_nonequi));

        //abs_diff_lt_precision(demag.h(state_full), nonequi_demag.h(state_nonequi), 2e-3);
        //abs_diff_lt_precision(demag.h(state_full), nonequi_demag.h(state_nonequi), 3e-3);

        //rel_diff_lt_precision(demag.h(state_full), nonequi_demag.h(state_nonequi), 5e-4);
        //rel_diff_lt_precision(demag.h(state_full), nonequi_demag.h(state_nonequi), 6e-4);
        //rel_diff_lt_precision(demag.h(state_full), nonequi_demag.h(state_nonequi), 7e-4);

        //double reldiff = rel_diff_upperbound(demag.h(state_full), nonequi_demag.h(state_nonequi), false, 1e0, 0.1, 0.9);
        //double absdiff = abs_diff_upperbound(demag.h(state_full), nonequi_demag.h(state_nonequi), false, 1e0, 0.1, 0.9);
        //std::cout << "rel diff = " << reldiff << std::endl;
        //std::cout << "abs diff = " << absdiff << std::endl;

        std::cout << "heff mean abs diff = " << mean_abs_diff(demag.h(state_full), nonequi_demag.h(state_nonequi)) << std::endl;
        std::cout << "heff mean rel diff = " << mean_rel_diff(demag.h(state_full), nonequi_demag.h(state_nonequi)) << std::endl << std::endl;

        std::cout << "N[000] abs diff = " << max_abs_diff(demag.todel_N(0, 0, 0), nonequi_demag.todel_N(0, 0, 0)) << std::endl;
        std::cout << "N[000] rel diff = " << max_rel_diff(demag.todel_N(0, 0, 0), nonequi_demag.todel_N(0, 0, 0)) << std::endl << std::endl;

        std::cout << "heff max abs diff = " << max_abs_diff(demag.h(state_full), nonequi_demag.h(state_nonequi)) << std::endl;
        std::cout << "heff max rel diff = " << max_rel_diff(demag.h(state_full), nonequi_demag.h(state_nonequi)) << std::endl << std::endl;

        std::cout << "N mean abs diff = " << mean_abs_diff(demag.todel_N, nonequi_demag.todel_N) << std::endl;
        std::cout << "N mean rel diff = " << mean_rel_diff(demag.todel_N, nonequi_demag.todel_N) << std::endl << std::endl;

        std::cout << "N max abs diff = " << max_abs_diff(demag.todel_N, nonequi_demag.todel_N) << std::endl;
        std::cout << "N max rel diff = " << max_rel_diff(demag.todel_N, nonequi_demag.todel_N) << std::endl;
    }
    return 0;
}
