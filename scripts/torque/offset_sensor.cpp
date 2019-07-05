#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace magnumaf;


int main(int argc, char** argv)
{
    // Checking input variables and setting GPU Device
    //timer total_time = af::timer::start();
    for (int i=0; i<argc; i++){std::cout << "Parameter " << i << " was " << argv[i] << std::endl;}
    std::string filepath(argc>1? argv[1]: "output_magnum.af/");
    setDevice(argc>2? std::stoi(argv[2]):0);
    info();

    // Parameter initialization
    const double x=1.e-9, y=1e-9, z=0.6e-9;
    const int nx = 2, ny=2 , nz=2;
    //const int nx = 1, ny=1 , nz=1;

    //Generating Objects
    Mesh mesh(nx, ny, nz, x/nx, y/ny, z/nz);
    Material material = Material();
    state.Ms    = 1.4e6;
    material.A     = 30e-12;
    material.alpha = 0.02;
    material.Ku1 = 1e4;
    material.Ku1_axis[0]=1;
    material.Ku1_axis[1]=0;
    material.Ku1_axis[2]=0;

    // Initial magnetic field
    af::array m = af::constant(0.0, nx, ny, nz, 3, f64);
    m(af::span, af::span, af::span, 0 )=1.;

    af::array pol = af::constant(0.0, nx, ny, nz, 3, f64);
    pol(af::span, af::span, af::span, 0 )=1/sqrt(2);
    pol(af::span, af::span, af::span, 1 )=1/sqrt(2);

    State state(mesh, material, m);
    vti_writer_micro(state.m, mesh , (filepath + "minit").c_str());

    LlgTerms llgterm;
    llgterm.push_back( LlgTerm (new DemagField(mesh, material, true)));
    llgterm.push_back( LlgTerm (new ExchangeField(mesh, material)));
    llgterm.push_back( LlgTerm (new SpinTransferTorqueField(pol, .3, .4, 2e10)));
    llgterm.push_back( LlgTerm (new UniaxialAnisotropyField(mesh, material)));
    LLGIntegrator Llg(llgterm);

    std::ofstream stream;
    stream.precision(12);
    stream.open ((filepath + "m.dat").c_str());

    // Relax
    timer t = af::timer::start();
    while (state.t < 1e-7){
        Llg.step(state);
        state.calc_mean_m(stream);
        //state.calc_mean_m(std::cout);
        //std::cout << afvalue(m(0, 0, 0, 0)) << " " << afvalue(m(0, 0, 0, 1)) << " " << afvalue(m(0, 0, 0, 2)) << " " << std::endl;
    }
    std::cout<<"timerelax [af-s]: "<< af::timer::stop(t) <<std::endl;
    vti_writer_micro(state.m, mesh , (filepath + "relax").c_str());
    return 0;
}
