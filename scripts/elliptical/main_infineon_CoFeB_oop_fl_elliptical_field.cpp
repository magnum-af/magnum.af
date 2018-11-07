#include "arrayfire.h"
#include "magnum_af.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>

double t_full_rotation = 800e-9;
double A = 0.025/(4e-7 * M_PI); // 0.025 T is half the linear range 
double B = 1.0 * A; // TODO
af::array zee_func(State state){
    double phi = 2 * M_PI * (state.t / t_full_rotation);
    array zee = constant(0.0,state.mesh.n0,state.mesh.n1,state.mesh.n2,3,f64);
    zee(span,span,span,0)=constant( A * std::cos(phi) ,state.mesh.n0,state.mesh.n1,state.mesh.n2,1,f64);
    zee(span,span,span,1)=constant( B * std::sin(phi) ,state.mesh.n0,state.mesh.n1,state.mesh.n2,1,f64);
    zee(span,span,span,2)=constant(-A * std::cos(phi) ,state.mesh.n0,state.mesh.n1,state.mesh.n2,1,f64);
    return  zee;
}

af::array zee_func_for_relax_in_init(State state){
    double phi = 0;
    array zee = constant(0.0,state.mesh.n0,state.mesh.n1,state.mesh.n2,3,f64);
    zee(span,span,span,0)=constant( A * std::cos(phi) ,state.mesh.n0,state.mesh.n1,state.mesh.n2,1,f64);
    zee(span,span,span,1)=constant( B * std::sin(phi) ,state.mesh.n0,state.mesh.n1,state.mesh.n2,1,f64);
    zee(span,span,span,2)=constant(-A * std::cos(phi) ,state.mesh.n0,state.mesh.n1,state.mesh.n2,1,f64);
    return  zee;
}
  
int main(int argc, char** argv)
{
    std::string filepath(argc >= 1? argv[1]: "../Data/Testing");
    if( argc >= 1 ){ filepath.append("/");}
    if( argc >= 2 ){ setDevice( std::stoi( argv[2]));}
    std::string path_mrelax(argc>3? argv[3]: "");
    std::cout<<"Writing into path "<<filepath.c_str()<<std::endl;
    std::cout.precision(24);
    info();


    // Parameter initialization
    const double x=800e-9, y=800e-9, z=1.3e-9;//[m] // Physical dimensions
    const int nx = 250, ny=250 ,nz=1;
    //const int nx = 400, ny=400 ,nz=1;
  
    //Generating Objects
    Mesh mesh(nx,ny,nz,x/nx,y/ny,z/nz);
    mesh.print(std::cout);
    Param param = Param();
    param.ms    = 1.58/param.mu0;//[J/T/m^3] == [Joule/Tesla/meter^3] = 1.75 T/mu_0
    param.A     = 15e-12;//[J/m]
    param.Ku1 = 1.3e-3/z;
    param.alpha = 0.02;

    long int n_cells=0;//Number of cells with Ms!=0
    State state(mesh,param, mesh.ellipse(n_cells, 2));

    vti_writer_micro(state.Ms, mesh ,(filepath + "Ms").c_str());
    vti_writer_micro(state.m, mesh ,(filepath + "minit").c_str());
    mesh.print(std::cout);

    std::vector<LlgTerm> llgterm;
    llgterm.push_back( LlgTerm (new DemagSolver(mesh,param)));
    llgterm.push_back( LlgTerm (new ExchSolver(mesh,param)));
    llgterm.push_back( LlgTerm (new ANISOTROPY(mesh,param)));
    llgterm.push_back( LlgTerm (new Zee(& zee_func_for_relax_in_init)));
    NewLlg Llg(llgterm);

    // Relaxation
    if(!exists (path_mrelax)){
        Llg.relax(state, 1e-7);
        vti_writer_micro(state.m, mesh ,(filepath + "mrelax").c_str());
        state.t=0; // Setting t=0 for hysteresis
    }
    else{
        std::cout << "found mrelax. loading magnetization" << std::endl;
        vti_reader(state.m, mesh, path_mrelax);
    }

    std::cout << "n_cells= " << n_cells << ", should be a*b*M_PI*mesh.n2= " << mesh.n0/2*mesh.n1/2*M_PI*mesh.n2 << std::endl;

    std::ofstream stream;
    stream.precision(12);
    stream.open ((filepath + "m.dat").c_str());
    stream << "# t	<mx>" << std::endl;
    state.calc_mean_m(stream, n_cells, Llg.llgterms[Llg.llgterms.size()-1]->h(state)(0,0,0,af::span));// To checkback H_zee for init
    Llg.llgterms.pop_back(); // Remove init zee field 

    timer t_hys = af::timer::start();
    Llg.llgterms.push_back( LlgTerm (new Zee(&zee_func))); //Rate in T/s
    while (state.t < t_full_rotation){
         Llg.step(state);
         state.calc_mean_m(stream, n_cells, Llg.llgterms[Llg.llgterms.size()-1]->h(state)(0,0,0,af::span));
         if( state.steps % 2000 == 0){
             vti_writer_micro(state.m, mesh ,(filepath + "m_hysteresis_"+std::to_string(state.steps)+std::to_string(state.t)).c_str());
         }
    }

    stream.close();
    std::cout<<"time full hysteresis [af-s]: "<< af::timer::stop(t_hys) <<std::endl;
    return 0;
}
