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
    //if(argc>1) setDevice(std::stoi(argv[2]));
    info();

    // Parameter initialization
    double length = 90e-9; //[nm]
    const double dx=0.25e-9;
    const int nx = (int)(length/dx);
    std::cout << "nx = "<< nx << std::endl;
  
    //Generating Objects
    Mesh mesh(nx,nx,1,dx,dx,dx);
    Param param = Param();
    param.ms    = 580000;
    param.A     = 15e-12;
    param.alpha = 1;
    param.D=3e-3;
    param.Ku1=0.6e6;
  
    param.J_atom=2.*param.A*dx;
    param.D_atom= param.D * pow(dx,2);
    param.K_atom=param.Ku1*pow(dx,3);
    param.p=param.ms*pow(dx,3);//Compensate nz=1 instead of nz=4
  
     // Initial magnetic field
    array m = constant(0.0,mesh.n0,mesh.n1,mesh.n2,3,f64);
    m(span,span,span,0)=1.;

    State state(mesh,param, m);
    //vti_writer_atom(state.m, mesh ,(filepath + "minit").c_str());
  
    std::vector<llgt_ptr> llgterm;
    llgterm.push_back( llgt_ptr (new ATOMISTIC_DEMAG(mesh)));
    llgterm.push_back( llgt_ptr (new ATOMISTIC_EXCHANGE(mesh)));
    llgterm.push_back( llgt_ptr (new ATOMISTIC_DMI(mesh,param)));
    llgterm.push_back( llgt_ptr (new ATOMISTIC_ANISOTROPY(mesh,param)));
    //llgterm.push_back( llgt_ptr (new DemagSolver(mesh,param)));
    //llgterm.push_back( llgt_ptr (new ExchSolver(mesh,param)));
    //llgterm.push_back( llgt_ptr (new DMI(mesh,param)));
    //llgterm.push_back( llgt_ptr (new ANISOTROPY(mesh,param)));
  
    LLG Llg(state,llgterm);

    const double E_inplane = Llg.E(state);
    state.m = constant(0.0,mesh.n0,mesh.n1,mesh.n2,3,f64);
    state.m(span,span,span,2)=1.;
    const double E_out_of_plane = Llg.E(state);
    const double Keff = pow(param.ms,2) * param.mu0 * pow(mesh.dx,3) * pow(mesh.n0,2);
    std::cout << E_inplane << "\t" << E_out_of_plane << "\t" << E_inplane - E_out_of_plane << "\t" << Keff << std::endl;
    std::cout << "dE/Keff = " << (E_inplane - E_out_of_plane)/Keff << std::endl;
    // TODO values?

    return 0;
}
