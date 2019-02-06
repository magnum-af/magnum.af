#include "arrayfire.h"
#include "magnum_af.hpp"

int main(int argc, char** argv)
{
    // Checking input variables and setting GPU Device
    timer total_time = af::timer::start();
    for (int i=0; i<argc; i++){cout << "Parameter " << i << " was " << argv[i] << std::endl;}
    info();
    
    // Parameter initialization
    const double x=5.e-7, y=1.25e-7, z=3.e-9;
    const int nx = 100, ny=25 ,nz=1;
    
    //Generating Objects
    Mesh mesh(nx,ny,nz,x/nx,y/ny,z/nz);
    Param param = Param();
    param.ms    = 8e5;
    param.A     = 1.3e-11;
    param.alpha = 1;
    
    // Initial magnetic field
    State state(mesh,param, mesh.init_sp4());
    
    LlgTerms llgterm;
    llgterm.push_back( LlgTerm (new ExchSolver(mesh,param)));
    NewLlg Llg(llgterm);
    
    af::array h = Llg.llgterms[0]->h(state);
    std::cout << "E= " << Llg.llgterms[0]->E(state) << std::endl;
    std::cout << "Eh=" << Llg.llgterms[0]->E(state, h) << std::endl;
    return 0;
}
