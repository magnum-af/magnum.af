#include "arrayfire.h"
#include "magnum_af.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>

double rate = 0.34e6 ; //[T/s]
double hzee_max = 0.12; //[T]

af::array zee_func(State state){
    double field_Tesla = 0;
    if(state.t < hzee_max/rate) field_Tesla = rate *state.t; 
    else if (state.t < 3*hzee_max/rate) field_Tesla = -rate *state.t + 2*hzee_max; 
    else if(state.t < 4*hzee_max/rate) field_Tesla = rate*state.t - 4*hzee_max; 
    else {field_Tesla = 0; std::cout << "WARNING ZEE time out of range" << std::endl;}
    array zee = constant(0.0,state.mesh.n0,state.mesh.n1,state.mesh.n2,3,f64);
    zee(span,span,span,0)=constant(field_Tesla/state.param.mu0 ,state.mesh.n0,state.mesh.n1,state.mesh.n2,1,f64);
    return  zee;
}
  
int main(int argc, char** argv)
{
    std::string filepath(argc >= 1? argv[1]: "../Data/Testing");
    if( argc >= 1 ){ filepath.append("/");}
    if( argc >= 2 ){ setDevice( std::stoi( argv[2]));}
    std::string path_mrelax(argc>3? argv[3]: "");
    std::cout<<"Writing into path "<<filepath.c_str()<<std::endl;
    std::cout.precision(24);
    info();

    // Parameter initialization
    const int nx = 250, ny=250 ,nz=1; // Discretization
    const double x=1600e-9, y=1600e-9, z=65e-9;//[m] // Physical dimensions
  
    //Generating Objects
    Mesh mesh(nx,ny,nz,x/nx,y/ny,z/nz);
    Param param = Param();
    param.ms    = 1.393e6;//[J/T/m^3] == [Joule/Tesla/meter^3] = 1.75 T/mu_0
    param.A     = 1.5e-11;//[J/m]
    param.alpha = 0.02;

    long int n_cells=0;//Number of cells with Ms!=0
    State state(mesh, param, mesh.init_vortex(n_cells));
    vti_writer_micro(state.Ms, mesh ,(filepath + "Ms").c_str());

    std::vector<LlgTerm> llgterm;
    llgterm.push_back( LlgTerm (new DemagSolver(mesh,param)));
    llgterm.push_back( LlgTerm (new ExchSolver(mesh,param)));
    NewLlg Llg(llgterm);

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
    state.calc_mean_m(stream, n_cells, Llg.llgterms[Llg.llgterms.size()-1]->h(state)(0,0,0,af::span));

    timer t_hys = af::timer::start();
    Llg.llgterms.push_back( LlgTerm (new Zee(&zee_func))); //Rate in T/s
    while (state.t < 4* hzee_max/rate){
         Llg.step(state);
         state.calc_mean_m(stream, n_cells, Llg.llgterms[Llg.llgterms.size()-1]->h(state)(0,0,0,af::span));
         if( state.steps % 2000 == 0){
             vti_writer_micro(state.m, mesh ,(filepath + "m_hysteresis_"+std::to_string(state.steps)).c_str());
         }
    }
    stream.close();
    std::cout<<"time full hysteresis [af-s]: "<< af::timer::stop(t_hys) <<std::endl;
    return 0;
}