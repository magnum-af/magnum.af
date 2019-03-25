#include "magnum_af.hpp"
#include "llg_terms/micro_demag.hpp"
#include "llg_terms/micro_nonequi_demag.hpp"

unsigned int zero_if_equal(af::array first, af::array second){
    unsigned int zero_if_equal = afvalue_u32(af::sum(af::sum(af::sum(af::sum(first != second, 0), 1), 2), 3));
        if (!zero_if_equal) std::cout << "\33[1;32mSucess:\33[0m zero_if_equal = " << zero_if_equal << std::endl;
        else std::cout << "\33[1;31mError!\33[0m zero_if_equal =" << zero_if_equal << std::endl;
    return zero_if_equal;
}


int main(int argc, char** argv)
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
    
    //zero_if_equal(demag.Nfft, nonequi_demag.Nfft);
    //zero_if_equal(demag.Nfft(af::span, af::span, 2, af::span), nonequi_demag.Nfft);
    return 0;
}
