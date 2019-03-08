#include "arrayfire.h"
#include "magnum_af.hpp"

int main(int argc, char** argv)
{
    // Checking input variables and setting GPU Device
    timer total_time = af::timer::start();
    for (int i=0; i<argc; i++){cout << "Parameter " << i << " was " << argv[i] << std::endl;}
    std::string filepath(argc>1? argv[1]: "output_magnum.af/");
    setDevice(argc>2? std::stoi(argv[2]):0);
    info();
    
    // Parameter initialization
    const double x=5.e-7, y=1.25e-7, z=3.e-9;
    const int nx = 100, ny=25 ,nz=1;
    
    //Generating Objects
    Mesh mesh(nx,ny,nz,x/nx,y/ny,z/nz);
    Material material = Material();
    material.ms    = 8e5;
    material.A     = 1.3e-11;
    material.alpha = 1;
    
    // Initial magnetic field
    State state(mesh,material, mesh.init_sp4());
    vti_writer_micro(state.m, mesh ,(filepath + "minit").c_str());
    
    LlgTerms llgterm;
    llgterm.push_back( LlgTerm (new DemagField(mesh,material, true)));
    llgterm.push_back( LlgTerm (new ExchangeField(mesh,material)));
    LLGIntegrator Llg(llgterm);
    
    std::ofstream stream;
    stream.precision(12);
    stream.open ((filepath + "m.dat").c_str());
    
    // Relax
    timer t = af::timer::start();
    while (state.t < 1e-9){
        Llg.step(state);
        state.calc_mean_m(stream);
    }
    std::cout<<"timerelax [af-s]: "<< af::timer::stop(t) <<std::endl;
    vti_writer_micro(state.m, mesh ,(filepath + "relax").c_str());

    // Prepare switch
    array zeeswitch = constant(0.0,1,1,1,3,f64);
    zeeswitch(0,0,0,0)=-24.6e-3/constants::mu0;
    zeeswitch(0,0,0,1)=+4.3e-3/constants::mu0;
    zeeswitch(0,0,0,2)=0.0;
    zeeswitch = tile(zeeswitch,mesh.n0,mesh.n1,mesh.n2);
    Llg.llgterms.push_back( LlgTerm (new Zee(zeeswitch)));
    state.material.alpha=0.02;

    // Switch
    t = af::timer::start();
    while (state.t < 2e-9){
        Llg.step(state);
        state.calc_mean_m(stream);
    }
    std::cout<<"time integrate 1ns [af-s]: "<< af::timer::stop(t) <<std::endl;
    vti_writer_micro(state.m, mesh ,(filepath + "2ns").c_str());
    stream.close();
    std::cout<<"total [af-s]: "<< af::timer::stop(total_time) <<std::endl;
    return 0;
}
