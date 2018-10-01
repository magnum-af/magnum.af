#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace af; 
typedef std::shared_ptr<LLGTerm> llgt_ptr; 

int main(int argc, char** argv)
{

    std::cout<<"argc = "<<argc<<std::endl;
     for (int i=0; i<argc; i++)
          cout << "Parameter " << i << " was " << argv[i] << "\n";
    
    std::string filepath(argc>1? argv[1]: "../Data/skyrmion_stoch");
    if(argc>0)filepath.append("/");
    std::cout<<"Writing into path "<<filepath.c_str()<<std::endl;

    setDevice(argc>2? std::stoi(argv[2]):0);
    info();

    // Parameter initialization
    const int nx = 90, nz=1;
    const double dx=1.0e-9;
    const double dz=0.6e-9;
  
    //Generating Objects
    Mesh mesh(nx,nx,nz,dx,dx,dz);
    Param param = Param();
    param.ms    = 580000;
    param.A     = 15e-12;
    param.alpha = 1;
    param.D=3e-3;
    param.Ku1=0.6e6;
  
     // Initial magnetic field
     array m = constant(0.0,mesh.n0,mesh.n1,mesh.n2,3,f64);
     m(span,span,span,2) = -1;
     for(int ix=0;ix<mesh.n0;ix++){
         for(int iy=0;iy<mesh.n1;iy++){
             const double rx=double(ix)-mesh.n0/2.;
             const double ry=double(iy)-mesh.n1/2.;
             const double r = sqrt(pow(rx,2)+pow(ry,2));
             if(r>nx/4.) m(ix,iy,span,2)=1.;
         }
     }
  
    State state(mesh,param, m);
    vti_writer_atom(state.m, mesh ,(filepath + "minit").c_str());
  
    std::vector<llgt_ptr> llgterm;
    //llgterm.push_back( llgt_ptr (new DemagSolver(mesh,param)));
    llgterm.push_back( llgt_ptr (new ExchSolver(mesh,param)));
    llgterm.push_back( llgt_ptr (new DMI(mesh,param)));
    llgterm.push_back( llgt_ptr (new ANISOTROPY(mesh,param)));
    
    
    LLG Llg(state,llgterm);

    std::cout << "mrelax.vti not found, starting relaxation" << std::endl;
    timer t = af::timer::start();
    while (state.t < 15.e-10){
        state.m=Llg.llgstep(state);
    }
    double timerelax= af::timer::stop(t);
    vti_writer_atom(state.m, mesh ,filepath + "relax");
  
    std::cout<<"timerelax [af-s]: "<< timerelax << " for "<<Llg.counter_accepted+Llg.counter_reject<<" steps, thereof "<< Llg.counter_accepted << " Steps accepted, "<< Llg.counter_reject<< " Steps rejected" << std::endl;
    return 0;
}
