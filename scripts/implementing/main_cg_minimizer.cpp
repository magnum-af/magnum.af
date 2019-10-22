#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace magnumafcpp;


int main(int argc, char** argv)
{
    std::string filepath(argc > 1 ? argv[1]: "../Data/Testing/");
    const double x=5.e-7, y=1.25e-7, z=3.e-9;
    const int nx = 100, ny=25 , nz=1;

    //Generating Objects
    Mesh mesh(nx, ny, nz, x/nx, y/ny, z/nz);
    Material material = Material();
    state.Ms    = 8e5;
    material.A     = 1.3e-11;

    // Initial magnetic field
    State state(mesh, material, mesh.init_sp4());
    vti_writer_micro(state.m, mesh , filepath + "minit");

    af::timer timer_llgterms = af::timer::start();
    CG_Minimizer minimizer = CG_Minimizer();
    minimizer.llgterms_.push_back( LlgTerm (new DemagField(mesh, material)));
    minimizer.llgterms_.push_back( LlgTerm (new ExchangeField(mesh, material)));
    std::cout<<"Llgterms assembled in "<< af::timer::stop(timer_llgterms) <<std::endl;

    timer_llgterms = af::timer::start();
    minimizer.Minimize(state);
    std::cout<<"Llgterms assembled in "<< af::timer::stop(timer_llgterms) <<std::endl;
    return 0;
}
