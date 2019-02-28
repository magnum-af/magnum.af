#include "arrayfire.h"
#include "magnum_af.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>

double rate = 0.20e6 ; //[T/s]
double hzee_max = 0.20; //[T]

af::array zee_func(State state){
    double field_Tesla = 0;
    if(state.t < hzee_max/rate) field_Tesla = rate *state.t; 
    else if (state.t < 3*hzee_max/rate) field_Tesla = -rate *state.t + 2*hzee_max; 
    else if(state.t < 4*hzee_max/rate) field_Tesla = rate*state.t - 4*hzee_max; 
    else {field_Tesla = 0; std::cout << "WARNING ZEE time out of range" << std::endl;}
    array zee = constant(0.0,state.mesh.n0,state.mesh.n1,state.mesh.n2,3,f64);
    zee(span,span,span,0)=constant(field_Tesla/state.material.mu0 ,state.mesh.n0,state.mesh.n1,state.mesh.n2,1,f64);
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
    const double x=800e-9, y=800e-9, z=1.3e-3/1.056e6;//[m] // z for 100mT lin range t_CoFeB = 1.3e-3/1.056e6  
    const int nx = 250, ny=250 ,nz=1;
  
    //Generating Objects
    Mesh mesh(nx,ny,nz,x/nx,y/ny,z/nz);
    mesh.print(std::cout);
    Material material = Material();
    material.ms    = 1.58/material.mu0;//[J/T/m^3] == [Joule/Tesla/meter^3] = 1.75 T/mu_0
    material.A     = 15e-12;//[J/m]
    material.Ku1   = 1.3e-3/z;// [J/m^3] // Ku1 = K_total - K_shape = Hk*Js/2/mu0 + Js^2/2/mu0 = | [Hk and Js in Tesla] | = ((0.1*1.58)/2/(4*pi*1e-7) + (1.58)^2/(2)/(4*pi*1e-7)) = 1.056e6
    material.alpha = 0.02;

    long int n_cells=0;//Number of cells with Ms!=0
    State state(mesh,material, mesh.ellipse(n_cells, 2));

    vti_writer_micro(state.Ms, mesh ,(filepath + "Ms").c_str());
    vti_writer_micro(state.m, mesh ,(filepath + "minit").c_str());
    mesh.print(std::cout);

    std::vector<LlgTerm> llgterm;
    llgterm.push_back( LlgTerm (new DemagField(mesh,material)));
    llgterm.push_back( LlgTerm (new ExchangeField(mesh,material)));
    llgterm.push_back( LlgTerm (new UniaxialAnisotropyField(mesh,material)));
    LLGIntegrator Llg(llgterm);

    // Relaxation
    if(!exists (path_mrelax)){
        Llg.relax(state, 1e-7);
        vti_writer_micro(state.m, mesh ,(filepath + "mrelax").c_str());
        state.t=0; // Setting t=0 for hysteresis
    }
    else{
        std::cout << "found mrelax. loading magnetization" << std::endl;
        vti_reader(state.m, mesh, path_mrelax);
    }

    std::cout << "n_cells= " << n_cells << ", should be a*b*M_PI*mesh.n2= " << mesh.n0/2*mesh.n1/2*M_PI*mesh.n2 << std::endl;

    std::ofstream stream;
    stream.precision(12);
    stream.open ((filepath + "m.dat").c_str());
    stream << "# t	<mx>" << std::endl;

    timer t_hys = af::timer::start();
    Llg.llgterms.push_back( LlgTerm (new Zee(&zee_func))); //Rate in T/s
    while (state.t < 4* hzee_max/rate){
         Llg.step(state);
         state.calc_mean_m(stream, n_cells, Llg.llgterms[Llg.llgterms.size()-1]->h(state)(0,0,0,af::span));
         if( state.steps % 2000 == 0){
             vti_writer_micro(state.m, mesh ,(filepath + "m_hysteresis_"+std::to_string(state.steps)));
         }
    }

    stream.close();
    std::cout<<"time full hysteresis [af-s]: "<< af::timer::stop(t_hys) <<std::endl;
    return 0;
}
