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
    const int nx = 100, ny=25 ,nz=1;
    
    Mesh mesh(nx,ny,nz,x/nx,y/ny,z/nz);
    Material material = Material();
    material.ms    = 8e5;
    material.A     = 1.3e-11;
    material.alpha = 1;

    DemagField demag = DemagField(mesh,material, true, false, 4);
    NonEquiDemagField nonequi_demag = NonEquiDemagField(mesh,material, true, false, 4);
    unsigned int zero_if_equal_var = zero_if_equal(demag.Nfft, nonequi_demag.Nfft);
    return 0;
}
