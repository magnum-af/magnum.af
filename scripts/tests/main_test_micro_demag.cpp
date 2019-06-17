#include "arrayfire.h"
#include "magnum_af.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <string>
#include <memory>
using namespace af; typedef std::shared_ptr<LLGTerm> llgt_ptr;
void calcm(State state, std::ostream& myfile);

bool compare(double a, double b){
    std::cout << "COM:"<< a <<"," << b <<","<<fabs(a-b)/fabs(a+b)<<std::endl;
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
    Material material = Material();
    state.Ms    = 1e5;

    // Initial magnetic field
    array m = constant(0.0,mesh.n0,mesh.n1,mesh.n2,3,f64);
    m(span,span,span,0) = 1;
    State state(mesh,material, m);

    //MICROMAGNETIC
    std::vector<llgt_ptr> llgterm;
    llgterm.push_back( llgt_ptr (new DemagField(mesh,material)));
    LLG Llg(state,llgterm);

    //std::cout.precision(12);
    //std::cout << "E_d_micro  = " << Llg.E(state) << std::endl;
    //std::cout << "Analytical = " << 1./6. * x * y * z * pow(state.Ms,2) * constants::mu0 << std::endl;

    if (compare(Llg.E(state),1./6. * x * y * z * pow(state.Ms,2) * constants::mu0)){
        std::cout << "!!! Test FAILED !!!" << std::endl;
        return 1;
    }
    //TODO  NOTE: Test fails for opencl (Error: 2.37e-8 > threshold), passes for cpu (Error: 4.05e-13 < threshold)
    return 0;
}
