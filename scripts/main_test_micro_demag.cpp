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

bool compare(double a, double b){
    //std::cout << "COM:"<< a <<"," << b <<","<<fabs(a-b)/fabs(a+b)<<std::endl;
    const double threshold = 1e-12;
    if(a == 0 && b == 0) return false;
    if(fabs(a-b)/fabs(a+b)< threshold) return false;
    else return true;
}

int main(int argc, char** argv)
{
    // Parameter initialization
    const double x=1.e-9, y=1.e-9, z=1.e-9;
    const int nx = 10, ny=10 ,nz=10;
  
    
    //Generating Objects
    Mesh mesh(nx,ny,nz,x/nx,y/ny,z/nz);
    Param param = Param();
    param.ms    = 1e5;
  
    // Initial magnetic field
    array m = constant(0.0,mesh.n0,mesh.n1,mesh.n2,3,f64);
    m(span,span,span,0) = 1;
    State state(mesh,param, m);
  
    //MICROMAGNETIC
    std::vector<llgt_ptr> llgterm;
    llgterm.push_back( llgt_ptr (new DemagSolver(mesh,param)));
    LLG Llg(state,llgterm);
  
    //std::cout.precision(12);
    //std::cout << "E_d_micro  = " << Llg.E(state) << std::endl;
    //std::cout << "Analytical = " << 1./6. * x * y * z * pow(param.ms,2) * param.mu0 << std::endl;
  
    if (compare(Llg.E(state),1./6. * x * y * z * pow(param.ms,2) * param.mu0)){
        std::cout << "!!! Test FAILED !!!" << std::endl;
        return 1;
    }
  
    return 0;
}
