#include "arrayfire.h"
#include "magnum_af.hpp"

int main(int argc, char** argv)
{
    std::string filepath(argc > 1 ? argv[1]: "../Data/Testing/");

    const double x=5.e-7, y=1.25e-7, z=3.e-9;
    const int nx = 100, ny=25 ,nz=1;
    
    //Generating Objects
    Mesh mesh(nx,ny,nz,x/nx,y/ny,z/nz);
    Param param = Param();
    param.ms    = 8e5;
    param.A     = 1.3e-11;
    param.alpha = 1;
    
    // Initial magnetic field
    State state(mesh, param, mesh.init_sp4());
    vti_writer_micro(state.m, mesh ,(filepath + "minit"));
    
    af::timer timer_llgterms = af::timer::start();
    LBFGS_Minimizer minimizer = LBFGS_Minimizer();
    minimizer.llgterms_.push_back( LlgTerm (new DemagSolver(mesh,param)));
    minimizer.llgterms_.push_back( LlgTerm (new ExchSolver(mesh,param)));
    std::cout<<"Llgterms assembled in [s]: "<< af::timer::stop(timer_llgterms) <<std::endl;

    double f = minimizer.Minimize(state);
    std::cout << "main: f= "<< f << std::endl;
    std::cout << green("Starting second run ")<< std::endl;
    f = minimizer.Minimize(state);
    std::cout << "main: f= "<< f << std::endl;
    af::print("minimizer", af::mean(state.m,0));
    vti_writer_micro(state.m, state.mesh, filepath+"m_minimized");
    NewLlg llg = NewLlg(minimizer.llgterms_);
    std::cout << red("E= ") << llg.E(state)  << green(" (as reference)") << std::endl;
    return 0;
}
