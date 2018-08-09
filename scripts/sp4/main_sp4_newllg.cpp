#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace af; typedef std::shared_ptr<LLGTerm> llgt_ptr; 

void calcm(State state, std::ostream& myfile);

int main(int argc, char** argv)
{
     std::cout<<"argc"<<argc<<std::endl;
     for (int i=0; i<argc; i++)
          cout << "Parameter " << i << " was " << argv[i] << "\n";
    
    std::string filepath(argc>1? argv[1]: "../../Data/Testing");
    if(argc>0)filepath.append("/");
    std::cout<<"Writing into path "<<filepath.c_str()<<std::endl;
    
    setDevice(argc>2? std::stoi(argv[2]):0);
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
    param.afsync  = false;
    param.T  = 300;
    
    // Initial magnetic field
    array m = constant(0.0,mesh.n0,mesh.n1,mesh.n2,3,f64);
    m(seq(1,end-1),span,span,0) = constant(1.0,mesh.n0-2,mesh.n1,mesh.n2,1,f64);
    m(0,span,span,1 ) = constant(1.0,1,mesh.n1,mesh.n2,1,f64);
    m(-1,span,span,1) = constant(1.0,1,mesh.n1,mesh.n2,1,f64);
    State state(mesh,param, m);
    vti_writer_micro(state.m, mesh ,(filepath + "minit").c_str());
    
    NewLlg Llg = NewLlg();
    Llg.llgterms.push_back( LlgTerm (new DemagSolver(mesh,param)));
    Llg.llgterms.push_back( LlgTerm (new ExchSolver(mesh,param)));
    
    std::ofstream stream;
    stream.precision(12);
    stream.open ((filepath + "m.dat").c_str());
    
    timer t = af::timer::start();
    while (state.t < 5.e-10){
        Llg.step(state);
        calcm(state,stream);
    }
    double timerelax= af::timer::stop(t);
    std::cout<<"timerelax [af-s]: "<< timerelax <<std::endl;
    vti_writer_micro(state.m, mesh ,(filepath + "relax").c_str());

    // Prepare switch
    array zeeswitch = constant(0.0,1,1,1,3,f64);
    zeeswitch(0,0,0,0)=-24.6e-3/param.mu0;
    zeeswitch(0,0,0,1)=+4.3e-3/param.mu0;
    zeeswitch(0,0,0,2)=0.0;
    zeeswitch = tile(zeeswitch,mesh.n0,mesh.n1,mesh.n2);
    Llg.llgterms.push_back( LlgTerm (new Zee(zeeswitch)));
    state.param.alpha=0.02;

    while (state.t < 1.5e-9){
        Llg.step(state);
        calcm(state,stream);
    }
    vti_writer_micro(state.m, mesh ,(filepath + "2ns").c_str());
    stream.close();
    return 0;
}

void calcm(State state, std::ostream& myfile){
    myfile << std::setw(12) << state.t << "\t" <<meani(state.m,0)<< "\t" <<meani(state.m,1)<< "\t" <<meani(state.m,2)<< "\t" << std::endl;
}
