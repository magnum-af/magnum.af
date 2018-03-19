#include "arrayfire.h"
#include "pth_mag.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <string>
#include <memory>
using namespace af; typedef std::shared_ptr<LLGTerm> llgt_ptr; 
void calcm(State state, std::ostream& myfile);
int main(int argc, char** argv)
{
    std::cout<<"argc"<<argc<<std::endl;
    for (int i=0; i<argc; i++) cout << "Parameter " << i << " was " << argv[i] << "\n";
    std::string filepath(argc>1? argv[1]: "../Data/Testing");
    if(argc>0)filepath.append("/");
    std::cout<<"Writing into path "<<filepath.c_str()<<std::endl;
    setDevice(argc>2? std::stoi(argv[2]):0);
    info();
  
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
    array m = constant(0.0,mesh.n0,mesh.n1,mesh.n2,3,f64);
    return 0;
}
