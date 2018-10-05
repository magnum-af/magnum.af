#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace af; 

int main(int argc, char** argv)
{
    std::string filepath(argc >= 1? argv[1]: "data");
    if( argc >= 1 ){ filepath.append("/");}
    if( argc >= 2 ){ setDevice( std::stoi( argv[2]));}
    std::cout<<"Writing into path "<<filepath.c_str()<<std::endl;
    std::string path_mrelax(argc>3? argv[3]: "");
    info();

    // Parameter initialization
    const int nxy = 90, nz=1;
    const double dx=1.0e-9;
    const double dz=0.6e-9;
  
    //Generating Objects
    Mesh mesh(nxy,nxy,nz,dx,dx,dz);
    Param param = Param();
    param.ms    = 580000;
    param.A     = 15e-12;
    param.alpha = 1;
    param.D=3e-3;
    param.Ku1=0.6e6;
  
    State state(mesh, param, mesh.skyrmconf());
    vti_writer_micro(state.m, mesh ,(filepath + "minit").c_str());
  
    NewLlg Llg;
    Llg.llgterms.push_back( LlgTerm (new DemagSolver(mesh,param)));
    Llg.llgterms.push_back( LlgTerm (new ExchSolver(mesh,param)));
    Llg.llgterms.push_back( LlgTerm (new DMI(mesh,param)));
    Llg.llgterms.push_back( LlgTerm (new ANISOTROPY(mesh,param)));
  
    if(!exists (path_mrelax)){
        std::cout << "mrelax.vti not found, starting relaxation" << std::endl;
        Llg.relax(state);
        vti_writer_micro(state.m, mesh ,filepath + "relax");
        state.t=0;
    }
    else{
        std::cout << "found mrelax. loading magnetization" << std::endl;
        vti_reader(state.m, state.mesh, path_mrelax);
    }

    array last   = constant( 0,mesh.dims,f64);
    last(span,span,span,2)=1;
    
    std::vector<State> inputimages; 
    inputimages.push_back(state);
    inputimages.push_back(State(mesh,param, last));
  
    String string(state, inputimages, 60, 5e-14, Llg.llgterms);
    string.run(filepath, 1e-13, 1e-28, 10000);
    return 0;
}
