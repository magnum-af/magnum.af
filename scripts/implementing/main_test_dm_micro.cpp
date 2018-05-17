#include "arrayfire.h"
#include "magnum_af.hpp"
using namespace af; 
typedef std::shared_ptr<LLGTerm> llgt_ptr; 
int main(int argc, char** argv)
{
    //const double x = 1.e-0, y = 1.e-0, z = 1.e-0;
    const int nx = 4, ny = 4 ,nz = 4;
    Mesh mesh(nx,ny,nz,1,1,1);
    Param param = Param();
    param.ms    = 1;
    param.D     = 1;
    param.D_axis[0] = 1;
    param.D_axis[1] = 0;
    param.D_axis[2] = 0;
    param.A     = 1;

    //array m = constant(0.0,mesh.n0,mesh.n1,mesh.n2,3,f64);
    //m(span,span,span,0) = constant(1.0,mesh.n0,mesh.n1,mesh.n2,1,f64);
    std::cout << mesh.n0 << std::endl;
    array m = constant (0, dim4(mesh.n0,mesh.n1,mesh.n2,3),f64);
    m(span,span,span,0)=1.;
    //array m = iota(dim4(mesh.n0,mesh.n1,mesh.n2,3),dim4(1,1,1,1),f64);
    print("m",m);
    //m=reorder(m,2,1,0,3);
    State state(mesh,param, m);
    std::vector<llgt_ptr> llgterm;
    llgterm.push_back( llgt_ptr (new DMI(mesh,param)));
    LLG Llg(state,llgterm);
    
    print("DMI", Llg.Fieldterms[0]->h(state));
    return 0;
}
