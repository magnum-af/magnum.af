#include "magnum_af.hpp"
#include "func.hpp"
#include "llg_terms/micro_demag_pthread.hpp"

int main(int argc, char** argv)
{
    for (int nz = 1; nz < 11; nz++){
        const double x=5.e-7, y=1.25e-7, z=3.e-9;
        const int nx = 400, ny=250;
        std::cout << "nz = " << nz << std::endl;
        //const int nx = 400, ny=250 ,nz=10;
        Mesh mesh(nx,ny,nz,x/nx,y/ny,z/nz);
        Material material = Material();
        DemagFieldMultithread PthreadDemag (mesh, material, true, false, 16);
        DemagField Demag (mesh, material, true, false);

        unsigned int zero_if_equal = afvalue_u32(af::sum(af::sum(af::sum(af::sum(Demag.Nfft != PthreadDemag.Nfft, 0), 1), 2), 3));
        if (!zero_if_equal) std::cout << "\33[1;32mSucess:\33[0m zero_if_equal = " << zero_if_equal << std::endl;
        else std::cout << "\33[1;31mError!\33[0m zero_if_equal =" << zero_if_equal << std::endl;
    }
    return 0;
}
    //af::print("Demag", Demag.Nfft);
    //af::print("Pthread", PthreadDemag.Nfft);
    //af::print("bool", af::sum(af::sum(af::sum(af::sum(boolean, 0), 1), 2), 3));
    //af::array boolean = ( Demag.Nfft != PthreadDemag.Nfft);
    //std::cout << boolean.type() << std::endl;
    //af::array zero_if_equal = af::sum(af::sum(af::sum(af::sum(boolean, 0), 1), 2), 3);
    //std::cout << zero_if_equal.type() << std::endl;
    //double zero_if_equal = afvalue(af::sum(af::sum(af::sum(af::sum(boolean, 0), 1), 2), 3));

