#include "../src/magnum_af.hpp"
#include "arrayfire.h"

using namespace magnumaf;


int main(int argc, char** argv)
{
    // Checking input variables and setting GPU Device
    for (int i=0; i<argc; i++){std::cout << "Parameter " << i << " was " << argv[i] << std::endl;}
    std::string filepath(argc>1? argv[1]: "output_magnum.af/");
    af::setDevice(argc>2? std::stoi(argv[2]):0);
    af::info();

    // Parameter initialization
    const float x=5.e-7, y=1.25e-7, z=3.e-9;
    const int nx = 100, ny=25 , nz=1;
    const float A = 1.3e-11;

    //Generating Objects
    Mesh mesh(nx, ny, nz, x/nx, y/ny, z/nz);

    // Initial magnetic field
    State state(mesh, 8e5, mesh.init_sp4());
    state.write_vti(filepath + "minit");

    auto demag = LlgTerm (new DemagField(mesh, true, true, 0));
    auto exch = LlgTerm (new ExchangeField(A));
    LLGIntegrator Llg(1, {demag, exch});

    std::ofstream stream;
    stream.precision(12);
    stream.open(filepath + "m.dat");

    // Relax
    StageTimer timer;
    while (state.t < 1e-9){
        Llg.step(state);
        state.calc_mean_m(stream);
    }
    timer.print_stage("relax ");
    state.write_vti(filepath + "relax");

    // Prepare switch
    af::array zeeswitch = af::constant(0.0, nx, ny, nz, 3, f32);
    zeeswitch(af::span, af::span, af::span, 0) = -24.6e-3/constants::mu0;
    zeeswitch(af::span, af::span, af::span, 1) = +4.3e-3/constants::mu0;
    Llg.llgterms.push_back( LlgTerm (new ExternalField(zeeswitch)));
    Llg.alpha = 0.02;

    // Switch
    while (state.t < 2e-9){
        Llg.step(state);
        state.calc_mean_m(stream);
    }
    state.write_vti(filepath + "2ns");
    stream.close();
    timer.print_stage("switch");
    timer.print_accumulated();
    return 0;
}
