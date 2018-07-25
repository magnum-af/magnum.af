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
    const int nx(argc>3? std::stoi(argv[3]):30);
    const int ny(argc>3? std::stoi(argv[3]):30);
    const int nz=1;
    const double dx(argc>4? std::stod(argv[4]):30);  
    double n_interp = 60;
    double string_dt=1e-13;
    const int string_steps = 10000;
  
  
    //Generating Objects
    Mesh mesh(nx,ny,nz,dx,dx,dx);
    Param param = Param();
    param.ms    = 1.1e6;
    param.A     = 1.6e-11;
    param.alpha = 1;
    param.D=2*5.76e-3;
    param.Ku1=6.4e6;
  
    param.J_atom=2.*param.A*dx;
    param.D_atom= param.D * pow(dx,2);
    param.K_atom=param.Ku1*pow(dx,3);
    param.p=param.ms*pow(dx,3);//Compensate nz=1 instead of nz=4
  
  
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
    
    //llgterm.push_back( llgt_ptr (new ATOMISTIC_DEMAG(mesh)));
    llgterm.push_back( llgt_ptr (new ATOMISTIC_EXCHANGE(mesh)));
    llgterm.push_back( llgt_ptr (new ATOMISTIC_DMI(mesh,param)));
    llgterm.push_back( llgt_ptr (new ATOMISTIC_ANISOTROPY(mesh,param)));
  
    LLG Llg(state,llgterm);
  
    timer t = af::timer::start();
    while (state.t < 8.e-10){
        state.m=Llg.llgstep(state);
    }
    double timerelax= af::timer::stop(t);
    vti_writer_atom(state.m, mesh ,(filepath + "relax").c_str());
  
    std::cout<<"timerelax [af-s]: "<< timerelax << " for "<<Llg.counter_accepted+Llg.counter_reject<<" steps, thereof "<< Llg.counter_accepted << " Steps accepted, "<< Llg.counter_reject<< " Steps rejected" << std::endl;
    return 0;
}
