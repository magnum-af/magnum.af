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
    const int steps_full_rotation =(argc > 5 ? std::stoi(argv[5]) : 200); 
    std::cout<<"Writing into path "<<filepath.c_str()<<std::endl;
    std::cout.precision(24);
    info();

    // Defining lamdas
    auto zee_func = [ steps_full_rotation, A, B ] ( State state ) -> af::array {
        double phi = 2. * M_PI * (double)state.steps / (double)steps_full_rotation;
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
    material.ms    = 1.75/constants::mu0;//[J/T/m^3] == [Joule/Tesla/meter^3] = 1.75 T/mu_0
    material.A     = 1.5e-11;//[J/m]
    std::cout << "A=" << A << "B= " << B << "steps_full_rotation=" << steps_full_rotation << std::endl;
    State state(mesh, material, mesh.init_vortex());
    vti_writer_micro(state.Ms, mesh ,(filepath + "Ms").c_str());

    Material param_stress = material;
    param_stress.Ku1 = 1400;
    param_stress.Ku1_axis[0]=1;
    param_stress.Ku1_axis[1]=0;
    param_stress.Ku1_axis[2]=0;
    param_stress.ms = material.ms;//TODO should be taken form state in the future

    af::timer timer_llgterms = af::timer::start();
    LBFGS_Minimizer minimizer(1e-6, 1000, 0);
    //LBFGS_Minimizer minimizer = LBFGS_Minimizer(1e-6, 1000, 0);// This fails on GTO with current gcc version
    minimizer.of_convergence.open(filepath + "minimizer_convergence.dat");
    minimizer.llgterms_.push_back( LlgTerm (new DemagField(mesh,material)));
    minimizer.llgterms_.push_back( LlgTerm (new ExchangeField(mesh,material)));
    minimizer.llgterms_.push_back( LlgTerm (new UniaxialAnisotropyField(mesh,param_stress)));
    minimizer.llgterms_.push_back( LlgTerm (new ExternalField(zee_func)));
    std::cout<<"Llgterms assembled in "<< af::timer::stop(timer_llgterms) <<std::endl;

    // Starting Hysteresis loop
    std::ofstream stream;
    stream.precision(18);
    stream.open ((filepath + "m.dat").c_str());
    stream << "# t	<mx>    <my>    <mz>    <h_x>    <h_y>    <h_z>" << std::endl;

    timer t_rot = af::timer::start();
    for (int i = 0; i <= steps_full_rotation; i++){
        minimizer.Minimize(state);
        state.calc_mean_m_steps(stream, minimizer.llgterms_.end()[-1]->h(state)(0,0,0,af::span));
        if( state.steps % 2000 == 0){
            vti_writer_micro(state.m, mesh ,filepath + "m_rotation_"+std::to_string(state.steps));
        }
        state.steps++;
    }
    stream.close();
    std::cout<<"time full rotation [af-s]: "<< af::timer::stop(t_rot) <<std::endl;
    return 0;
}
