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
    //double length = 00.1e-9; //[nm]
    const int nx = 2500;
    const double dx=1.e-9;
    //const int nx = (int)(length/dx);
    //std::cout << "nx = "<< nx << std::endl;
  
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
  
    std::vector<llgt_ptr> llgterm;
    llgterm.push_back( llgt_ptr (new ATOMISTIC_DEMAG(mesh)));
  
    LLG Llg(state,llgterm);

    // Inplane 
    vti_writer_micro(state.m, state.mesh ,(filepath + "m_inplane").c_str());
    const double E_inplane = Llg.E(state);
    Llg.write_fieldterms_micro(state, (filepath + "demagfield_in_plane").c_str());

    // Out-of-plane
    state.m = constant(0.0,mesh.n0,mesh.n1,mesh.n2,3,f64);
    state.m(span,span,span,2)=1.;
    vti_writer_micro(state.m, state.mesh ,(filepath + "m_ouf_of_plane").c_str());
    const double E_out_of_plane = Llg.E(state);
    Llg.write_fieldterms_micro(state, (filepath + "demagfield_out_of_plane").c_str());

    // Keff 
    const double Keff = - 0.5 * pow(param.ms,2) * param.mu0 * pow(mesh.dx,3) * mesh.n0 * mesh.n1 * mesh.n2;

    const double g_electron = 2.002319;
    const double correction = 1.07831;

    const double mub = 9.27400999e-24;
    const double deltaEd = 2 * M_PI * correction * pow( g_electron * mub,2)/pow(dx,3);
    //TODO
    //const double deltaEd = 2 * M_PI * correction * pow( g_electron * mub,2)/pow(dx,3) * mesh.n0 * mesh.n1 * mesh.n2;

    std::cout << "deltaEd = " << deltaEd << std::endl;

    std::cout << "E_inplance= " << E_inplane << "\n" << "E_out_of_plane= " << E_out_of_plane << "\n" << "E_inplane - E_out_of_plane= " << E_inplane - E_out_of_plane << "\n" << "Keff= " <<  Keff << std::endl;
    std::cout << "dE/Keff = " << (E_inplane - E_out_of_plane)/Keff << std::endl;

    return 0;
}
