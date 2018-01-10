//Following paper Journal of Magnetism and Magnetic Materials 233 (2001) 296–304 from Scholz, Schrefl, Fidler
//ATTENTION, ERROR IN PAPER EQ 22 second expression, should be: (epx(chi)-1)/(sqrt(pi*chi)*erfi(sqrt(chi)))
#include <iostream>
#include <complex>
#include "arrayfire.h"
#include "pth_mag.hpp"
using namespace af;
using namespace std::complex_literals;
using Faddeeva::erfi;
 typedef std::shared_ptr<LLGTerm> llgt_ptr; 
void calcm(State state, std::ostream& myfile);

//double pth_erfi(double x){
//    std::complex<double> cx = x ;
//    std::complex<double> im = 0.0 + 1i;
//    std::complex<double> result = - 1i * erf(1i*x);
//    return result.real;
//}

double mean_mz_analytical (double chi){
    return (exp(chi) - 1.)/(sqrt(M_PI) * sqrt(chi) * erfi(sqrt(chi))); //TODO erfi 
}

int main(int argc, char** argv)
{
    std::cout<<"argc"<<argc<<std::endl;
    for (int i=0; i<argc; i++){cout << "Parameter " << i << " was " << argv[i] << "\n";}
    std::string filepath(argc>1? argv[1]: "../Data/Testing");
    if(argc>0)filepath.append("/");
    std::cout<<"Writing into path "<<filepath.c_str()<<std::endl;
  
    setDevice(argc>2? std::stoi(argv[2]):0);
    info();
  
    // Parameter initialization
    const double x=1.e-9, y=1.e-9, z=1.e-9;
    const int nx = 1, ny=1 ,nz=1;
    double dt = 1e-14;
  
    //Generating Objects
    Mesh mesh(nx,ny,nz,x/nx,y/ny,z/nz);
    Param param = Param();
    param.gamma = 2.211e5;
    param.ms    = 1281197;
    param.Ku1   = 6.9e6;
    param.alpha = 0.1;
    param.T = 10;

    //Integration param
    //unsigned long int relax_steps   = 1e2;
    unsigned long int measure_steps = 2e4;
  
    // Initial magnetic field
    array m = constant(0.0,mesh.n0,mesh.n1,mesh.n2,3,f64);
    m(0,0,0,2)=1.;
    State state(mesh,param, m);
    vti_writer_micro(state.m, mesh ,(filepath + "minit").c_str());
  
    std::vector<llgt_ptr> llgterm;
    llgterm.push_back( llgt_ptr (new ANISOTROPY(mesh,param)));
    Stochastic_LLG Stoch(state,llgterm,dt,"Heun");
  
    std::ofstream stream;
    stream.precision(12);
    stream.open ((filepath + "m.dat").c_str());
    
    ////RELAX
    //timer t = af::timer::start();
    //for (unsigned long int i = 0; i < relax_steps; i++){
    //    Stoch.step(state,dt); 
    //    calcm(state,stream);
    //}
    //double timerelax= af::timer::stop(t);
    //std::cout<<"timerelax [af-s]: "<< timerelax <<std::endl;
    //vti_writer_micro(state.m, mesh ,(filepath + "relax").c_str());


    //timer t2 = af::timer::start();
    //double timemeasure= af::timer::stop(t2);
    //std::cout<<"timemeasure [af-s]: "<< timemeasure <<std::endl;

    //MEASURE
    Stoch.param.T  = 10;
    double mean_mz=0;
    for (unsigned long int i = 0; i < measure_steps; i++){
        Stoch.step(state,dt); 
        mean_mz+=afvalue(state.m(0,0,0,2));
        //if(i%1000==0) std::cout<< i <<"  "<< mean_mz <<std::endl;
        calcm(state,stream);
    }

    mean_mz=mean_mz/measure_steps;
    double chi = (state.param.Ku1 * state.mesh.V) / (state.param.kb * Stoch.param.T);
    std::cout<<"Analytical <mz>/Ms  = "<< mean_mz_analytical(chi) << " at "<<Stoch.param.T<<" K" <<std::endl;
    std::cout<<"Calculated <mz>/Ms  = "<< mean_mz <<"\n"<<std::endl;

    Stoch.param.T  = 50;
    //state.param.T  = 50;
    mean_mz=0;
    for (unsigned long int i = 0; i < measure_steps; i++){
        Stoch.step(state,dt); 
        mean_mz+=afvalue(state.m(0,0,0,2));
        //if(i%1000==0) std::cout<< i <<"  "<< mean_mz <<std::endl;
        calcm(state,stream);
    }

    mean_mz=mean_mz/measure_steps;
    chi = (state.param.Ku1 * state.mesh.V) / (state.param.kb * Stoch.param.T);
    std::cout<<"Analytical <mz>/Ms  = "<< mean_mz_analytical(chi) << " at "<<Stoch.param.T<<" K" <<std::endl;
    std::cout<<"Calculated <mz>/Ms  = "<< mean_mz <<"\n"<<std::endl;

    Stoch.param.T  = 200;
    //state.param.T  = 200;
    mean_mz=0;
    for (unsigned long int i = 0; i < measure_steps; i++){
        Stoch.step(state,dt); 
        mean_mz+=afvalue(state.m(0,0,0,2));
        //if(i%1000==0) std::cout<< i <<"  "<< mean_mz <<std::endl;
        calcm(state,stream);
    }

    mean_mz=mean_mz/measure_steps;
    chi = (state.param.Ku1 * state.mesh.V) / (state.param.kb * Stoch.param.T);
    std::cout<<"Analytical <mz>/Ms  = "<< mean_mz_analytical(chi) << " at "<<Stoch.param.T<<" K" <<std::endl;
    std::cout<<"Calculated <mz>/Ms  = "<< mean_mz <<"\n"<<std::endl;

    //std::cout<<"Reference T=10        = "<< 0.98979 <<std::endl;
    //std::cout<<"Reference T=50        = "<< 0.94268 <<std::endl;
    //std::cout<<"Reference T=200       = "<< 0.71976 <<std::endl;
  
    stream.close();
    return 0;
}
void calcm(State state, std::ostream& myfile){
    myfile << std::setw(12) << state.t << "\t" <<meani(state.m,0)<< "\t" <<meani(state.m,1)<< "\t" <<meani(state.m,2)<< "\t" << std::endl;
}
