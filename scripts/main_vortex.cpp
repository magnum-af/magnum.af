#include "arrayfire.h"
#include "pth_mag.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <string>
#include <memory>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
using namespace af; typedef std::shared_ptr<LLGTerm> llgt_ptr; 
void calcm(State state, std::ostream& myfile);
inline bool exists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}
  
int main(int argc, char** argv)
{
    std::cout<<"argc"<<argc<<std::endl;
    for (int i=0; i<argc; i++) cout << "Parameter " << i << " was " << argv[i] << "\n";
    std::string filepath(argc>1? argv[1]: "../Data/Testing");
    if(argc>0)filepath.append("/");
    std::cout<<"Writing into path "<<filepath.c_str()<<std::endl;
    setDevice(argc>2? std::stoi(argv[2]):0);
    info();
    std::cout.precision(12);

    // Parameter initialization
    const int nx = 250, ny=250 ,nz=1;
    const double x=1600e-9, y=1600e-9, z=65e-9;//[m] // Physical dimensions
  
    //Generating Objects
    Mesh mesh(nx,ny,nz,x/nx,y/ny,z/nz);
    Param param = Param();
    param.ms    = 1.393e6;//[J/T/m^3] == [Joule/Tesla/meter^3] = 1.75 T/mu_0
    param.A     = 1.5e-11;//[J/m]
    param.alpha = 0.02;
    // Initial magnetic field
    array Ms = constant(0.0,mesh.n0,mesh.n1,mesh.n2,3,f64);
    for(int ix=0;ix<mesh.n0;ix++){
        for(int iy=0;iy<mesh.n1;iy++){
            const double rx=double(ix)-mesh.n0/2.;
            const double ry=double(iy)-mesh.n1/2.;
            const double r = sqrt(pow(rx,2)+pow(ry,2));
            if(r<nx/2.){
                Ms(ix,iy,span,span)=param.ms;
            }
        }
    }


    array m;
    State state(mesh,param, m);
    state.Ms=Ms;
    vti_writer_micro(state.Ms, mesh ,(filepath + "Ms").c_str());

    std::vector<llgt_ptr> llgterm;
    llgterm.push_back( llgt_ptr (new DemagSolver(mesh,param)));
    llgterm.push_back( llgt_ptr (new ExchSolver(mesh,param)));
    //llgterm.push_back( llgt_ptr (new Zee(zee,mesh,param)));
    LLG Llg(state,llgterm);

    //obtaining relaxed magnetization
    if(!exists (filepath+"mrelax.vti")){
        std::cout << "mrelax.vti not found, starting relaxation" << std::endl;
        m = constant(0.0,mesh.n0,mesh.n1,mesh.n2,3,f64);
        for(int ix=0;ix<mesh.n0;ix++){
            for(int iy=0;iy<mesh.n1;iy++){
                const double rx=double(ix)-mesh.n0/2.;
                const double ry=double(iy)-mesh.n1/2.;
                const double r = sqrt(pow(rx,2)+pow(ry,2));
                if(r<nx/2.){
                    Ms(ix,iy,span,span)=param.ms;
                    if(r==0.){
                        m(ix,iy,span,2)= 1;
                    }
                    else{
                        m(ix,iy,span,0)=-ry/r;
                        m(ix,iy,span,1)= rx/r;
                        m(ix,iy,span,2)= sqrt(nx)/r;
                    }
                }
            }
        }
        vti_writer_micro(m, mesh ,(filepath + "minit_nonnormalized").c_str());
        m=renormalize_handle_zero_values(m);
        state.m=m;
        vti_writer_micro(state.m, mesh ,(filepath + "minit_renorm").c_str());
        timer t = af::timer::start();
        double E_prev=1e20;
        double E=0;
        while (fabs((E_prev-E)/E_prev) > 1e-11){
            state.m=Llg.llgstep(state);
            if( state.steps % 100 == 0){
                E_prev=E;
                E=Llg.E(state);
                std::cout << "step " << state.steps << " time= " << state.t << " E= " << E << " fabs= " << fabs((E_prev-E)/E_prev) << std::endl;
                vti_writer_micro(state.m, mesh ,(filepath + "m_current_relaxation"+std::to_string(state.steps)).c_str());
            }
        }
        std::cout<<"timerelax [af-s]: "<< af::timer::stop(t) <<std::endl;
        vti_writer_micro(state.m, mesh ,(filepath + "mrelax").c_str());
        std::cout << "state.t= "<< state.t << std::endl;
        state.t=0;//setting new zero time
    }
    else{
        std::cout << "found mrelax. loading magnetization" << std::endl;
        vti_reader(m, mesh, filepath+"mrelax.vti");
        state.m=m;
    }
    vti_writer_micro(state.m, mesh ,(filepath + "todel_check_mrelax").c_str());


    //Calc n_cells:
    long int n_cells=0;//Number of cells with Ms!=0
    for(int ix=0;ix<mesh.n0;ix++){
        for(int iy=0;iy<mesh.n1;iy++){
            const double rx=double(ix)-mesh.n0/2.;
            const double ry=double(iy)-mesh.n1/2.;
            const double r = sqrt(pow(rx,2)+pow(ry,2));
            if(r<nx/2.){
                n_cells++;
            }
        }
    }
    std::cout << "n_cells" << n_cells << std::endl;
    //array zee = constant(0.0,1,1,1,3,f64);
    //zee(0,0,0,0)=120e-3/param.mu0;//120 mT field
    //zee = tile(zee,mesh.n0,mesh.n1,mesh.n2);


    //while (state.t < 1.e-8){
    //    state.m=Llg.llgstep(state);
    //    if( state.steps % 100 == 0){
    //        std::cout << "step " << state.steps << " time= " << state.t << " E= " << Llg.E(state) <<  std::endl;
    //        vti_writer_micro(state.m, mesh ,(filepath + "m_current_relaxation"+std::to_string(state.steps)).c_str());
    //    }
    //}

  
    //Rate of Zeeman change: 700/1500 0.34 mT/ns
    // Simulated 1400 ns

    return 0;
}
