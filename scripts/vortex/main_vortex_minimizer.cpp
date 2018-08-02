#include "arrayfire.h"
#include "magnum_af.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <string>
#include <memory>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
using namespace af; typedef std::shared_ptr<LLGTerm> llgt_ptr; 

void calc_mean_m(const State& state, const long int n_cells,  std::ostream& myfile){
    array sum_dim3 = sum(sum(sum(state.m,0),1),2);
    myfile << std::setw(12) << state.t << "\t" << afvalue(sum_dim3(span,span,span,0))/n_cells << "\t" << afvalue(sum_dim3(span,span,span,1))/n_cells<< "\t" << afvalue(sum_dim3(span,span,span,2))/n_cells << std::endl;
}

void calc_mean_m(const State& state, const long int n_cells,  std::ostream& myfile, double hzee){
    array sum_dim3 = sum(sum(sum(state.m,0),1),2);
    myfile << std::setw(12) << state.t << "\t" << afvalue(sum_dim3(span,span,span,0))/n_cells << "\t" << afvalue(sum_dim3(span,span,span,1))/n_cells<< "\t" << afvalue(sum_dim3(span,span,span,2))/n_cells << "\t" << hzee << std::endl;
}

double hzee_max = 0.12; //[T]
int quater_steps=100; // One 4th of total steps

af::array zee_func(State state){
    double field_Tesla = 0;
    double rate = hzee_max/quater_steps; //[T/s]
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
    std::cout<<"argc"<<argc<<std::endl;
    for (int i=0; i<argc; i++) cout << "Parameter " << i << " was " << argv[i] << "\n";
    std::string filepath(argc>1? argv[1]: "../Data/Testing");
    if(argc>0)filepath.append("/");
    std::cout<<"Writing into path "<<filepath.c_str()<<std::endl;
    setDevice(argc>2? std::stoi(argv[2]):0);
    std::string path_mrelax(argc>3? argv[3]: "");
    info();
    std::cout.precision(24);

    // Parameter initialization
    const int nx = 250, ny=250 ,nz=1;
    //const int nx = 400, ny=400 ,nz=1;
    const double x=1600e-9, y=1600e-9, z=65e-9;//[m] // Physical dimensions
  
    //Generating Objects
    Mesh mesh(nx,ny,nz,x/nx,y/ny,z/nz);
    Param param = Param();
    param.ms    = 1.393e6;//[J/T/m^3] == [Joule/Tesla/meter^3] = 1.75 T/mu_0
    param.A     = 1.5e-11;//[J/m]
    param.alpha = 0.02;
    // Initial magnetic field
    array Ms = constant(0.0,mesh.n0,mesh.n1,mesh.n2,3,f64);
    array m  = constant(0.0,mesh.n0,mesh.n1,mesh.n2,3,f64);
    for(int ix=0;ix<mesh.n0;ix++){
        for(int iy=0;iy<mesh.n1;iy++){
            const double rx=double(ix)-mesh.n0/2.;
            const double ry=double(iy)-mesh.n1/2.;
            const double r = sqrt(pow(rx,2)+pow(ry,2));
            if(r<nx/2.){
                Ms(ix,iy,span,span)=param.ms;
                if(r==0.){
                    m(ix,iy,span,2)= 1;
                }
                else{
                    m(ix,iy,span,0)=-ry/r;
                    m(ix,iy,span,1)= rx/r;
                    m(ix,iy,span,2)= sqrt(nx)/r;
                }
            }
        }
    }
    State state(mesh,param, m);
    state.Ms=Ms;
    vti_writer_micro(state.Ms, mesh ,(filepath + "Ms").c_str());

    vti_writer_micro(state.m, mesh ,(filepath + "minit_nonnormalized").c_str());
    state.m=renormalize_handle_zero_values(state.m);
    vti_writer_micro(state.m, mesh ,(filepath + "minit_renorm").c_str());

    af::timer timer_llgterms = af::timer::start();
    Minimizer minimizer("BB", 1e-10, 1e-5, 1e4, 100);
    minimizer.llgterms.push_back( LlgTerm (new DemagSolver(mesh,param)));
    minimizer.llgterms.push_back( LlgTerm (new ExchSolver(mesh,param)));
    minimizer.llgterms.push_back( LlgTerm (new ANISOTROPY(mesh,param)));
    std::cout<<"Llgterms assembled in "<< af::timer::stop(timer_llgterms) <<std::endl;

    //obtaining relaxed magnetization
    if(!exists (path_mrelax)){
        timer t = af::timer::start();
        std::cout<<"timerelax [af-s]: "<< af::timer::stop(t) <<std::endl;
        minimizer.minimize(state);
        vti_writer_micro(state.m, mesh ,(filepath + "mrelax").c_str());
    }
    else{
        std::cout << "found mrelax. loading magnetization" << std::endl;
        vti_reader(state.m, state.mesh, filepath+"mrelax.vti");
    }
    vti_writer_micro(state.m, mesh ,(filepath + "todel_check_mrelax").c_str());

    //Calc n_cells:
    long int n_cells=0;//Number of cells with Ms!=0
    for(int ix=0;ix<mesh.n0;ix++){
        for(int iy=0;iy<mesh.n1;iy++){
            const double rx=double(ix)-mesh.n0/2.;
            const double ry=double(iy)-mesh.n1/2.;
            const double r = sqrt(pow(rx,2)+pow(ry,2));
            if(r<nx/2.){
                n_cells++;
            }
        }
    }
    std::cout << "n_cells= " << n_cells << ", should be nx^2*M_PI/4.= " << pow(nx,2)*M_PI/4. << std::endl;
    calc_mean_m(state,n_cells,std::cout);

    std::ofstream stream;
    stream.precision(12);
    stream.open ((filepath + "m.dat").c_str());
    stream << "# t	<mx>    <my>    <mz>    hzee" << std::endl;
    calc_mean_m(state,n_cells,stream);

    timer t_hys = af::timer::start();
    double rate = hzee_max/quater_steps; //[T/s]
    minimizer.llgterms.push_back( LlgTerm (new Zee(&zee_func)));
    while (state.t < 4* hzee_max/rate){
        minimizer.minimize(state);
        calc_mean_m(state,n_cells,stream,afvalue(minimizer.llgterms[3]->h(state)(0,0,0,0)));
        state.t+=1.;
        state.steps++;
        if( state.steps % 10 == 0){
            vti_writer_micro(state.m, mesh ,(filepath + "m_hysteresis_"+std::to_string(state.steps)).c_str());
        }
    }
    stream.close();
    std::cout<<"time full hysteresis [af-s]: "<< af::timer::stop(t_hys) <<std::endl;
    return 0;
}
