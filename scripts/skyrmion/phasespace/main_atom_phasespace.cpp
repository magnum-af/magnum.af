#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace af; typedef std::shared_ptr<LLGTerm> llgt_ptr; 
int main(int argc, char** argv)
{
    std::string filepath(argc >= 1? argv[1]: "data");
    if( argc >= 1 ){ filepath.append("/");}
    if( argc >= 2 ){ setDevice( std::stoi( argv[2]));}
    std::cout<<"Writing into path "<<filepath.c_str()<<std::endl;
    std::string path_mrelax(argc>3? argv[3]: "");
    info();
    
    // Parameter initialization
    const int nx = 112, ny=112 ,nz=1;//nz=5 -> lz=(5-1)*dx
    const double dx=2.715e-10;
  
  
    double n_interp = 60;
    double string_dt=5e-14;
    const int string_max_steps = 10000;
    double rel_diff = 1e-12;
    double abs_diff = 1e-27;

  
    //Generating Objects
    Mesh mesh(nx,ny,nz,dx,dx,dx);
    Param param = Param();
    param.ms    = 1.1e6;
    param.A     = 1.6e-11;
    param.alpha = 1;
    param.afsync  = false;
    param.D = (argc >= 3 ? std::stod(argv[3]) : 0.01152);
    param.Ku1=(argc >= 4 ? std::stod(argv[4]) : 6400000);

    param.set_atomistic_from_micromagnetic(mesh.dx);
  
    std::cout<<"D="<<param.D<<std::endl;
    std::cout<<"Ku1="<<param.Ku1<<std::endl;
    std::cout<<"D_atom="<<param.D_atom<<std::endl;
    std::cout<<"Ku1_atom="<<param.K_atom<<std::endl;
  
    State state(mesh, param, mesh.skyrmconf());
    vti_writer_micro(state.m, mesh ,(filepath + "minit").c_str());
  
    NewLlg Llg;
    Llg.llgterms.push_back( LlgTerm (new ATOMISTIC_DEMAG(mesh)));
    Llg.llgterms.push_back( LlgTerm (new ATOMISTIC_EXCHANGE(mesh)));
    Llg.llgterms.push_back( LlgTerm (new ATOMISTIC_DMI(mesh,param)));
    Llg.llgterms.push_back( LlgTerm (new ATOMISTIC_ANISOTROPY(mesh,param)));

    Llg.relax(state);
    vti_writer_micro(state.m, mesh ,filepath + "relax");
    state.t=0;
  
  
    array last   = constant( 0,mesh.dims,f64);
    last(span,span,span,2)=1;

    std::vector<State> inputimages; 
    inputimages.push_back(state);
    inputimages.push_back(State(mesh,param, last));
  
    String string(state, inputimages, n_interp, string_dt , Llg.llgterms);
    string.run(filepath, rel_diff, abs_diff, string_max_steps);
    return 0;
}
