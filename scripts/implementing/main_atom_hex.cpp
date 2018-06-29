#include "magnum_af.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <string>
#include <memory>
using namespace af; typedef std::shared_ptr<LLGTerm> llgt_ptr; 
int main(int argc, char** argv)
{
    std::cout<<"argc"<<argc<<std::endl;
     for (int i=0; i<argc; i++)
          cout << "Parameter " << i << " was " << argv[i] << "\n";
    
    std::string filepath(argc>1? argv[1]: "../Data/Testing");
    if(argc>0)filepath.append("/");
    std::cout<<"Writing into path "<<filepath.c_str()<<std::endl;

    setDevice(argc>2? std::stoi(argv[2]):0);
    info();

    // Parameter initialization
    const int nx = 5, ny=5 ,nz=1;//nz=5 -> lz=(5-1)*dx
    const double dx=1.e-10;
    //Generating Objects
    Mesh mesh(nx,ny,nz,dx,dx,dx);
    Param param = Param();
    param.p    = 1.;
    //param.alpha = 1.;
    param.J_atom=1;

     // Initial magnetic field
     array m = constant(0.0,mesh.n0,mesh.n1,mesh.n2,3,f64);
     //m(span,span,span,2) = -1;
     m(1,1,0) = 1;
  
    State state(mesh,param, m);
    vti_writer_atom(state.m, mesh ,(filepath + "minit").c_str());
  
    std::vector<llgt_ptr> llgterm;
    
    llgterm.push_back( llgt_ptr (new ATOMISTIC_EXCHANGE(mesh)));
    af::print("llgterm",llgterm[0]->h(state));
    //LLG Llg(state,llgterm);
  
    return 0;
}
