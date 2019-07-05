#include "arrayfire.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <string>
#include <memory>
#include "magnum_af.hpp"

using namespace magnumaf;


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
    const double x=1, y=1, z=1;
    //const double x=5.e-7, y=1.25e-7, z=3.e-9;
    const int nx = 10, ny=10 , nz=1;

    //Generating Objects
    Mesh mesh(nx, ny, nz, x/nx, y/ny, z/nz);
    Material material = Material();
    material.p    = 1;
    material.alpha = 1;
    material.J_atom= 1;

    // Initial magnetic field
    array m = constant(0.0, mesh.n0, mesh.n1, mesh.n2, 3, f64);
    State state(mesh, material, m);
    vti_writer_atom(state.m, mesh , (filepath + "minit").c_str());

    std::vector<llgt_ptr> llgterm;
    llgterm.push_back( llgt_ptr (new AtomisticExchangeField(mesh)));
    std::cout << "test " << std::endl;
    LLG Llg(state, llgterm);
    std::cout << "test " << std::endl;
    //TODO: not calling step or llgterm->h causes segfault in cleanup
    //af::print("h", llgterm[0]->h(state));
    //state.m=Llg.step(state);
    return 0;
}
