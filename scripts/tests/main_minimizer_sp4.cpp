#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace af;
int main(int argc, char** argv)
{
    std::string filepath(argc>1? argv[1]: "../../Data/Testing");
    if(argc>0)filepath.append("/");
    setDevice(argc>2? std::stoi(argv[2]):0);
    std::string path_mrelax(argc>3? argv[3]: "");
    // Parameter initialization
    const double x=5.e-7, y=1.25e-7, z=3.e-9;
    const int nx = 100, ny=25 ,nz=1;
    
    //Generating Objects
    Mesh mesh(nx,ny,nz,x/nx,y/ny,z/nz);
    Param param = Param();
    param.ms    = 8e5;
    param.A     = 1.3e-11;
    param.alpha = 1;
    
    // Initial magnetic field
    array m = constant(0.0,mesh.n0,mesh.n1,mesh.n2,3,f64);
    m(seq(1,end-1),span,span,0) = constant(1.0,mesh.n0-2,mesh.n1,mesh.n2,1,f64);
    m(0,span,span,1 ) = constant(1.0,1,mesh.n1,mesh.n2,1,f64);
    m(-1,span,span,1) = constant(1.0,1,mesh.n1,mesh.n2,1,f64);
    State state(mesh,param, m);
    if(exists (path_mrelax)){
        std::cout << "found mrelax. loading magnetization" << std::endl;
        vti_reader(state.m, mesh, path_mrelax);
    }

    LBFGS_Minimizer minimizer = LBFGS_Minimizer(1e-6, 1000, 0);
    minimizer.llgterms_.push_back( LlgTerm (new DemagSolver(mesh,param)));
    minimizer.llgterms_.push_back( LlgTerm (new ExchSolver(mesh,param)));
    std::cout << "E= " << minimizer.llgterms_[1]->E(state) << std::endl;
    std::cout << "Eh=" << minimizer.llgterms_[1]->E(state, minimizer.llgterms_[1]->h(state)) << std::endl;
        
    State state2(mesh,param, m);
    LBFGS_Minimizer minimizer2 = LBFGS_Minimizer(1e-6, 1000, 0);
    minimizer2.llgterms_.push_back( LlgTerm (new DemagSolver(mesh,param)));
    minimizer2.llgterms_.push_back( LlgTerm (new ExchSolverSparse(mesh,param)));
    std::cout << "E= " << minimizer2.llgterms_[1]->E(state) << std::endl;
    std::cout << "Eh=" << minimizer2.llgterms_[1]->E(state, minimizer2.llgterms_[1]->h(state)) << std::endl;
    return 0;
}
