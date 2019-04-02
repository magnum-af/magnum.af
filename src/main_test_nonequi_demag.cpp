#include "magnum_af.hpp"
#include "llg_terms/micro_demag.hpp"
#include "llg_terms/micro_nonequi_demag.hpp"

/////< Checks absolute value of point-wise difference: returns 0 if | x - y | < precision, 1 otherwise
//unsigned int zero_if_equal(af::array first, af::array second, double precision = 4e-12, bool verbose = true){
//    unsigned int zero_if_equal = afvalue_u32(af::sum(af::sum(af::sum(af::sum( !(af::abs(first - second) < precision), 0), 1), 2), 3));
//    if (verbose){
//        if (!zero_if_equal) std::cout << "\33[1;32mSucess:\33[0m zero_if_equal = " << zero_if_equal << std::endl;
//        else std::cout << "\33[1;31mError!\33[0m zero_if_equal =" << zero_if_equal << std::endl;
//    }
//    return zero_if_equal;
//}
    //unsigned int zero_if_equal = afvalue_u32(af::sum(af::sum(af::sum(af::sum(first != second, 0), 1), 2), 3));
    //af::array temp = af::abs(first - second);
    //af::print("temp", temp);
    //af::print("temp", !(temp < precision));
    //af::print("zero_if_equal", first);
    //af::print("zero_if_equal", second);
    //af::print("zero_if_equal", first != second);
    //if(first.isreal()) std::cout << "isreal" << std::endl;
    //if(first.isbool()) std::cout << "isbool" << std::endl;
    //if(first.isrealfloating()) std::cout << "isrealfloating" << std::endl;
    //std::cout.precision(24);
    //std::cout << "first=" << afvalue(first(0, 0, 0, 1)) << ", second=" << afvalue(second(0, 0, 0, 1))<< std::endl;


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

        af::print("full h", demag.h(state_full));

        af::print("nonequi_demag H", nonequi_demag.h(state_nonequi));

        abs_diff_lt_precision(demag.h(state_full), nonequi_demag.h(state_nonequi), 8e-4);
        abs_diff_lt_precision(demag.h(state_full), nonequi_demag.h(state_nonequi), 4e-3);

        rel_diff_lt_precision(demag.h(state_full), nonequi_demag.h(state_nonequi), 5e-4);
        rel_diff_lt_precision(demag.h(state_full), nonequi_demag.h(state_nonequi), 6e-4);
        rel_diff_lt_precision(demag.h(state_full), nonequi_demag.h(state_nonequi), 7e-4);

        double prec = 1e0;
        double prec_prev = prec;
        while(rel_diff_lt_precision(demag.h(state_full), nonequi_demag.h(state_nonequi), prec, false))
        {
            std::cout << "prec = " <<prec << std::endl;
            prec_prev = prec;
            prec = 0.1 * prec;
        }
        //std::cout << "finished width prec_prev = " <<prec_prev << std::endl;
        //std::cout << "finished width prec = " <<prec << std::endl;

        prec = prec_prev;
        while(rel_diff_lt_precision(demag.h(state_full), nonequi_demag.h(state_nonequi), prec, false))
        {
            std::cout << "prec = " <<prec << std::endl;
            prec_prev = prec;
            prec = 0.9 * prec;
        }
        std::cout << "finished width prec_prev = " <<prec_prev << " First false with " << prec << std::endl;
        rel_diff_lt_precision(demag.h(state_full), nonequi_demag.h(state_nonequi), prec_prev, true);
    }
    return 0;
}
