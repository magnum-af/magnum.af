//Following paper Journal of Magnetism and Magnetic Materials 233 (2001) 296–304 from Scholz, Schrefl, Fidler
//ATTENTION, ERROR IN PAPER EQ 22 second expression, should be: (epx(chi)-1)/(sqrt(pi*chi)*erfi(sqrt(chi)))
//Integrals in EQ 22 only consider positive z, so we take fabs(mz) ! (this compensates switching which would lead to an average around zero)
#include <iostream>
#include <complex>
#include "arrayfire.h"
#include "pth_mag.hpp"
using namespace af;
using namespace std::complex_literals;
using Faddeeva::erfi;
typedef std::shared_ptr<LLGTerm> llgt_ptr;

//Mathematica:
//(e^x-1)/(sqrt(pi*x)*erfi(sqrt(x))) =(int(exp(x * z^2)*z) dz from 0 to 1 )/(int(exp(x * z^2)) dz from 0 to 1)
double mean_mz_analytical (double chi){
    return (exp(chi) - 1.)/(sqrt(M_PI) * sqrt(chi) * erfi(sqrt(chi))); //TODO erfi 
}

void set_m_to_z(State& state){
    state.m(span,span,span,0)=0.;
    state.m(span,span,span,1)=0.;
    state.m(span,span,span,2)=1.;
}
void calcm(State state, std::ostream& myfile){
    myfile << std::setw(12) << state.t << "\t" <<meani(state.m,0)<< "\t" <<meani(state.m,1)<< "\t" <<meani(state.m,2)<< "\t" << std::endl;
}

int main(int argc, char** argv)
{
    std::cout<<"argc"<<argc<<std::endl;
    for (int i=0; i<argc; i++){cout << "Parameter " << i << " was " << argv[i] << "\n";}
    std::string filepath(argc>1? argv[1]: "../Data/rigid");
    if(argc>0)filepath.append("/");
    std::cout<<"Writing into path "<<filepath.c_str()<<std::endl;
    std::ofstream stream;
    stream.precision(12);
    stream.open ((filepath + "m.dat").c_str());
    std::ofstream stream2;
    stream2.precision(12);
    stream2.open ((filepath + "rigid.dat").c_str());
  
    setDevice(argc>2? std::stoi(argv[2]):0);
    info();
  
    //Integration param
    unsigned long int measure_steps = 6e4;

    stream2<<"#measure_steps = "<<measure_steps<<std::endl;
    stream2<<"# dt << Stoch.param.T<<  mean_mz <<  abs_mean_mz << mean_mz_analytical(chi) <<Stoch.param.T<<";
    stream2<<" mean_mz << abs_mean_mz << mean_mz_analytical(chi) <<Stoch.param.T<<  mean_mz";
    stream2<<" abs_mean_mz << mean_mz_analytical(chi)<<measure_steps"<<std::endl;
    // Parameter initialization
    const double x=1.e-9, y=1.e-9, z=1.e-9;
    const int nx = 1, ny=1 ,nz=1;
  
    //Generating Objects
    Mesh mesh(nx,ny,nz,x/nx,y/ny,z/nz);
    Param param = Param();
    param.gamma = 2.211e5;
    param.ms    = 1281197;
    param.Ku1   = 6.9e6;
    param.alpha = 0.1;
    param.T = 10;
  
    // Initial magnetic field
    array m = constant(0.0,mesh.n0,mesh.n1,mesh.n2,3,f64);
    m(0,0,0,2)=1.;
    State state(mesh,param, m);//ATTENTION, to be set in each loop
    std::vector<llgt_ptr> llgterm;
    Stochastic_LLG Stoch(state,llgterm,0.,"Heun");//ATTENTION, to be set in each loop
  
    //Declare Variables
    double mean_mz;
    double abs_mean_mz;
    double chi;
    //Initialize others
    double dt = 1e-15;

    //MEASURE
    
    for(int i=0;i<25;i++){
        std::cout<<"i= "<<i<<" dt = "<<dt<<std::endl;
        //T=10
        param.T=10.;
        state=State(mesh,param,m);
        llgterm.push_back( llgt_ptr (new ANISOTROPY(mesh,param)));
        Stoch = Stochastic_LLG(state,llgterm,dt,"Heun");
        llgterm.clear();
        mean_mz=0;
        abs_mean_mz=0;
        for (unsigned long int i = 0; i < measure_steps; i++){
            Stoch.step(state,dt); 
            mean_mz+=afvalue(state.m(0,0,0,2));
            abs_mean_mz+=fabs(afvalue(state.m(0,0,0,2)));
        }

        mean_mz=mean_mz/measure_steps;
        abs_mean_mz=abs_mean_mz/measure_steps;
        chi = (state.param.Ku1 * state.mesh.V) / (state.param.kb * Stoch.param.T);
        std::cout<<"at "<<Stoch.param.T<<" K" << std::endl; 
        std::cout<<"Calculated <mz>/Ms  = "<< mean_mz <<std::endl;
        std::cout<<"Cal  <fabs(mz)>/Ms  = "<< abs_mean_mz <<std::endl;
        std::cout<<"Analytical <mz>/Ms  = "<< mean_mz_analytical(chi) <<"\n" << std::endl;
        stream2<< std::setw(6)<< dt << "\t" << Stoch.param.T<< "\t" << mean_mz << "\t" << abs_mean_mz << "\t"<< mean_mz_analytical(chi)<< "\t";

        //T=50
        param.T=50.;
        state=State(mesh,param,m);
        llgterm.push_back( llgt_ptr (new ANISOTROPY(mesh,param)));
        Stoch = Stochastic_LLG(state,llgterm,dt,"Heun");
        llgterm.clear();
        mean_mz=0;
        abs_mean_mz=0;
        for (unsigned long int i = 0; i < measure_steps; i++){
            Stoch.step(state,dt); 
            mean_mz+=afvalue(state.m(0,0,0,2));
            abs_mean_mz+=fabs(afvalue(state.m(0,0,0,2)));
        }

        mean_mz=mean_mz/measure_steps;
        abs_mean_mz=abs_mean_mz/measure_steps;
        chi = (state.param.Ku1 * state.mesh.V) / (state.param.kb * Stoch.param.T);
        std::cout<<"at "<<Stoch.param.T<<" K" << std::endl; 
        std::cout<<"Calculated <mz>/Ms  = "<< mean_mz <<std::endl;
        std::cout<<"Cal  <fabs(mz)>/Ms  = "<< abs_mean_mz <<std::endl;
        std::cout<<"Analytical <mz>/Ms  = "<< mean_mz_analytical(chi) <<"\n" << std::endl;
        stream2<< Stoch.param.T << "\t"<< mean_mz << "\t" << abs_mean_mz << "\t"<< mean_mz_analytical(chi)<< "\t";

        //T=200
        param.T=200.;
        state=State(mesh,param,m);
        llgterm.push_back( llgt_ptr (new ANISOTROPY(mesh,param)));
        Stoch = Stochastic_LLG(state,llgterm,dt,"Heun");
        llgterm.clear();
        mean_mz=0;
        abs_mean_mz=0;
        for (unsigned long int i = 0; i < measure_steps; i++){
            Stoch.step(state,dt); 
            mean_mz+=afvalue(state.m(0,0,0,2));
            abs_mean_mz+=fabs(afvalue(state.m(0,0,0,2)));
        }

        mean_mz=mean_mz/measure_steps;
        abs_mean_mz=abs_mean_mz/measure_steps;
        chi = (state.param.Ku1 * state.mesh.V) / (state.param.kb * Stoch.param.T);
        std::cout<<"at "<<Stoch.param.T<<" K" << std::endl; 
        std::cout<<"Calculated <mz>/Ms  = "<< mean_mz <<std::endl;
        std::cout<<"Cal  <fabs(mz)>/Ms  = "<< abs_mean_mz <<std::endl;
        std::cout<<"Analytical <mz>/Ms  = "<< mean_mz_analytical(chi) <<"\n" << std::endl;
        stream2<< Stoch.param.T<< "\t" << mean_mz << "\t" << abs_mean_mz << "\t"<< mean_mz_analytical(chi)<<"\t"<<measure_steps<<std::endl;

        dt*=1.35;
    }

    stream.close();
    stream2.close();
    return 0;
}
