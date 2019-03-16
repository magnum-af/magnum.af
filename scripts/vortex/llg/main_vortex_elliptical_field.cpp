#include "arrayfire.h"
#include "magnum_af.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace af;

// argv[1]: filepath,  argv[2]: Device, argv[3]: A in mT, argv[4]: quotient B/A in percent, argv[5] t_full_rotation in seconds
int main(int argc, char** argv)
{
    std::string filepath(argc > 1? argv[1]: "../Data/Testing");
    if( argc > 1 ){ filepath.append("/");}
    setDevice( argc > 2 ? std::stoi( argv[2]) : 0);
    // Input a in mT, argv[3]=25 mT is converted to 0.025 T and divided by mu0
    const double A = double(argc > 3 ? std::stod(argv[3])*1e-3/(4e-7 * M_PI) : (double)(0.05/(4e-7 * M_PI)));
    const double B = double(argc > 4 ? std::stod(argv[4])/100 : 1.0) * A; // Input a in percent, B=1.0 == 100%
    const double t_full_rotation = double(argc > 5 ? std::stod(argv[5]) : (double)(800e-9));
    const std::string path_mrelax(argc>5? argv[5]: "");
    std::cout<<"Writing into path "<<filepath.c_str()<<std::endl;
    std::cout.precision(24);
    info();

    // Defining lamdas
    auto zee_func_for_relax_in_init= [ A, B ] ( State state ) -> af::array {
        double phi = 0;
        array zee = constant(0.0,state.mesh.n0,state.mesh.n1,state.mesh.n2,3,f64);
        zee(span,span,span,0)=constant( A * std::cos(phi) ,state.mesh.n0,state.mesh.n1,state.mesh.n2,1,f64);
        zee(span,span,span,1)=constant( B * std::sin(phi) ,state.mesh.n0,state.mesh.n1,state.mesh.n2,1,f64);
        zee(span,span,span,2)=constant( A * std::sin(phi) ,state.mesh.n0,state.mesh.n1,state.mesh.n2,1,f64);
        return  zee;
    };

    auto zee_func = [ t_full_rotation, A, B ] ( State state ) -> af::array {
        double phi = 2 * M_PI * (state.t / t_full_rotation);
        array zee = constant(0.0,state.mesh.n0,state.mesh.n1,state.mesh.n2,3,f64);
        zee(span,span,span,0)=constant( A * std::cos(phi) ,state.mesh.n0,state.mesh.n1,state.mesh.n2,1,f64);
        zee(span,span,span,1)=constant( B * std::sin(phi) ,state.mesh.n0,state.mesh.n1,state.mesh.n2,1,f64);
        zee(span,span,span,2)=constant( A * std::sin(phi) ,state.mesh.n0,state.mesh.n1,state.mesh.n2,1,f64);
        return  zee;
    };

    // Parameter initialization
    const int nx = 250, ny=250 ,nz=1; // Discretization
    const double x=1600e-9, y=1600e-9, z=65e-9;//[m] // Physical dimensions
  
    //Generating Objects
    Mesh mesh(nx,ny,nz,x/nx,y/ny,z/nz);
    Material material = Material();
    material.ms    = 1.393e6;//[J/T/m^3] == [Joule/Tesla/meter^3] = 1.75 T/mu_0
    material.A     = 1.5e-11;//[J/m]
    material.alpha = 0.02;

    std::cout << "A=" << A << "B= " << B << "t_full_rotation=" << t_full_rotation << std::endl;

    long int n_cells=0;//Number of cells with Ms!=0
    State state(mesh, material, mesh.init_vortex(n_cells));
    vti_writer_micro(state.Ms, mesh ,(filepath + "Ms").c_str());

    std::vector<LlgTerm> llgterm;
    llgterm.push_back( LlgTerm (new DemagField(mesh,material)));
    llgterm.push_back( LlgTerm (new ExchangeField(mesh,material)));
    llgterm.push_back( LlgTerm (new ExternalField( zee_func_for_relax_in_init)));
    LLGIntegrator Llg(llgterm);

    // Calculating relaxed initial magnetization or reading in given magnetization
    if(!exists (path_mrelax)){
        std::cout << "mrelax.vti not found, starting relaxation" << std::endl;
        vti_writer_micro(state.m, mesh ,(filepath + "minit_renorm").c_str());
        Llg.relax(state, 1e-7);
        vti_writer_micro(state.m, mesh ,(filepath + "mrelax").c_str());
        state.t=0; // Setting t=0 for hysteresis
    }
    else{
        std::cout << "found mrelax. loading magnetization" << std::endl;
        vti_reader(state.m, mesh, filepath+"mrelax.vti");
    }

    // Starting Hysteresis loop
    std::ofstream stream;
    stream.precision(12);
    stream.open ((filepath + "m.dat").c_str());
    stream << "# t	<mx>" << std::endl;
    state.calc_mean_m(stream, n_cells, Llg.llgterms[Llg.llgterms.size()-1]->h(state)(0,0,0,af::span));// To checkback H_zee for init
    Llg.llgterms.pop_back(); // Remove init zee field 

    timer t_hys = af::timer::start();
    Llg.llgterms.push_back( LlgTerm (new ExternalField(zee_func))); //Rate in T/s
    while (state.t < t_full_rotation){
         Llg.step(state);
         state.calc_mean_m(stream, n_cells, Llg.llgterms[Llg.llgterms.size()-1]->h(state)(0,0,0,af::span));
         if( state.steps % 2000 == 0){
             vti_writer_micro(state.m, mesh ,filepath + "m_hysteresis_"+std::to_string(state.steps));
         }
    }
    stream.close();
    std::cout<<"time full hysteresis [af-s]: "<< af::timer::stop(t_hys) <<std::endl;
    return 0;
}
