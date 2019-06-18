#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace af;
int main(int argc, char** argv)
{
    std::string filepath(argc>1? argv[1]: "../../Data/Testing");
    if(argc>0)filepath.append("/");
    // Parameter initialization
    const double x=5.e-7, y=1.25e-7, z=3.e-9;
    const int nx = 100, ny=25 , nz=1;

    //Generating Objects
    Mesh mesh(nx, ny, nz, x/nx, y/ny, z/nz);
    Material material = Material();
    state.Ms    = 8e5;
    material.A     = 1.3e-11;
    material.alpha = 1;

    // Initial magnetic field
    array m = constant(0.0, mesh.n0, mesh.n1, mesh.n2, 3, f64);
    m(seq(1, end-1), span, span, 0) = constant(1.0, mesh.n0-2, mesh.n1, mesh.n2, 1, f64);
    m(0, span, span, 1 ) = constant(1.0, 1, mesh.n1, mesh.n2, 1, f64);
    m(-1, span, span, 1) = constant(1.0, 1, mesh.n1, mesh.n2, 1, f64);
    State state(mesh, material, m);

    //LlgTerms llgterms;
    Minimizer minimizer;
    minimizer.llgterms.push_back( LlgTerm (new DemagField(mesh, material)));
    minimizer.llgterms.push_back( LlgTerm (new ExchangeField(mesh, material)));
    //minimizer.h(state);
    //af::print("h", minimizer.h(state));
    vti_writer_micro(state.m, mesh , (filepath + "init").c_str());
    minimizer.minimize(state);
    vti_writer_micro(state.m, mesh , (filepath + "minimized").c_str());
    std::cout<<"finished"<<std::endl;
    return 0;
}
