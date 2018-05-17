#include "arrayfire.h"
#include "magnum_af.hpp"
using namespace af; typedef std::shared_ptr<LLGTerm> llgt_ptr; 

bool compare(double a, double b){
    //std::cout << "COM:"<< a <<"," << b <<","<<fabs(a-b)/fabs(a+b)<<std::endl;
    if(a == 0 && b == 0) return false;
    if(fabs(a + b ) == 0) {std::cout << "DIVISION by 0 in compare" <<  std::endl; return true;}
    if(fabs(a-b)/fabs(a+b)<1e-15) return false;
    else return true;
}

int main(int argc, char** argv)
{
    info();
    std::cout.precision(32);
    int nx = 5, ny=5 ,nz=5;
    const double dx=2.715e-10;
  
    //Generating Objects
    Mesh mesh(nx,ny,nz,dx,dx,dx);
    Param param = Param();
    param.ms    = 1.1e6;
    param.A     = 1.6e-11;
    param.J_atom=2.*param.A*dx;
    param.p=param.ms*pow(dx,3);
  
    //-------------------------------------------------------
    array m = randu(mesh.n0,mesh.n1,mesh.n2,3,f64);
    State state(mesh,param, m);
  
    std::vector<llgt_ptr> llgterm;
    llgterm.push_back( llgt_ptr (new ATOMISTIC_EXCHANGE(mesh)));
    LLG Llg(state,llgterm);
    std::vector<llgt_ptr> llgterm2;
    llgterm2.push_back( llgt_ptr (new ExchSolver(mesh,param)));
    LLG Llg2(state,llgterm2);
    
    for (int x=0; x < nx; x++){
        for (int y=0; y < ny; y++){
            for (int z=0; z < nz; z++){
                for (int m=0; m <  3; m++){
                    if(compare(afvalue(Llg.Fieldterms[0]->h(state)(x,y,z,m)),afvalue(Llg2.Fieldterms[0]->h(state)(x,y,z,m)))){
                        std::cout <<"!!! TEST  FAILED at "<< x <<" " << y << " " << z << " " << m << std::endl;
                	std::cout << afvalue(Llg.Fieldterms[0]->h(state)(x,y,z,m)) << " , " << afvalue(Llg2.Fieldterms[0]->h(state)(x,y,z,m)) << std::endl;
                    }
                }
            }
        }
    }
    return 0;
}
