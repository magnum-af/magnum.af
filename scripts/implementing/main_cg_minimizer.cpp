#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace af;
int main()
{
    const double x=5.e-7, y=1.25e-7, z=3.e-9;
    const int nx = 100, ny=25 ,nz=1;
    
    //Generating Objects
    Mesh mesh(nx,ny,nz,x/nx,y/ny,z/nz);
    Param param = Param();
    param.ms    = 8e5;
    param.A     = 1.3e-11;
    param.alpha = 1;
    param.afsync  = false;
    param.T  = 300;
    
    // Initial magnetic field
    array m = constant(0.0,mesh.n0,mesh.n1,mesh.n2,3,f64);
    m(seq(1,end-1),span,span,0) = constant(1.0,mesh.n0-2,mesh.n1,mesh.n2,1,f64);
    m(0,span,span,1 ) = constant(1.0,1,mesh.n1,mesh.n2,1,f64);
    m(-1,span,span,1) = constant(1.0,1,mesh.n1,mesh.n2,1,f64);
    State state(mesh,param, m);
    //vti_writer_micro(state.m, mesh ,(filepath + "minit").c_str());
    
    af::timer timer_llgterms = af::timer::start();
    CG_Minimizer minimizer = CG_Minimizer();
    minimizer.llgterms.push_back( LlgTerm (new DemagSolver(mesh,param)));
    minimizer.llgterms.push_back( LlgTerm (new ExchSolver(mesh,param)));
    std::cout<<"Llgterms assembled in "<< af::timer::stop(timer_llgterms) <<std::endl;

    timer_llgterms = af::timer::start();
    minimizer.minimize(state);
    std::cout<<"Llgterms assembled in "<< af::timer::stop(timer_llgterms) <<std::endl;
    return 0;
}
