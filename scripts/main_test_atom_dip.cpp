#include "arrayfire.h"
#include "pth_mag.hpp"
using namespace af; typedef std::shared_ptr<LLGTerm> llgt_ptr; 

bool compare(double a, double b){
    if(fabs(a-b)/fabs(a+b)<1e-30) return false;
    else return true;
}

int main(int argc, char** argv)
{
    info();
    std::string filepath(argc>0? argv[1]: "../Data/Testing/");
    std::cout << filepath << std::endl;
    int nx = 2, ny=1 ,nz=1;//nz=5 -> lz=(5-1)*dx
    //const double dx=1;
    const double dx=2.715e-10;
  
    //Generating Objects
    Mesh mesh(nx,ny,nz,dx,dx,dx);
    Param param = Param();
    param.alpha = 1;
    param.afsync  = false;
    //param.p=1;
    param.p=9.274009994e-24;
  
    // Initial magnetic field
    //
    //-------------------------------------------------------
    array m = constant(0.0,mesh.n0,mesh.n1,mesh.n2,3,f64);
    m(0,0,0,0) = 0;
    m(0,0,0,1) = 0;
    m(0,0,0,2) = 1;

    m(1,0,0,0) = 0;
    m(1,0,0,1) = 0;
    m(1,0,0,2) = 1;
    State state(mesh,param, m);
    vti_writer_atom(state.m, mesh ,(filepath + "/minit").c_str());
  
    std::vector<llgt_ptr> llgterm;
    llgterm.push_back( llgt_ptr (new ATOMISTIC_DEMAG(mesh)));
    LLG Llg(state,llgterm);
    double analytical=- pow(param.p,2)*param.mu0/(4.*M_PI)/pow(dx,3);
    //std::cout << "ENERGY    = " << Llg.E(state) <<std::endl;
    std::cout << "Analytical= " << analytical <<std::endl;
    if(compare(Llg.E(state),analytical)) std::cout <<"!!! TEST FAILED !!!"<< std::endl;
    //std::cout << "Analytical= " << - pow(param.p,2)*param.mu0/(4.*M_PI)/pow(dx,3) <<std::endl;
    std::cout << "H_dip_1   = " <<0<<","<<0<<","<< param.p/(4*M_PI*pow(dx,3)) <<std::endl;
    std::cout << "H_dip_2   = " <<0<<","<<0<<","<< param.p/(4*M_PI*pow(dx,3)) <<std::endl;

    //-------------------------------------------------------
    m = constant(0.0,mesh.n0,mesh.n1,mesh.n2,3,f64);
    m(0,0,0,0) = 0;
    m(0,0,0,1) = 0;
    m(0,0,0,2) = 1;

    m(1,0,0,0) = 1;
    m(1,0,0,1) = 0;
    m(1,0,0,2) = 0;

    state.m=m;
    vti_writer_atom(state.m, mesh ,(filepath + "/minit2").c_str());
    std::cout << "ENERGY    = " << Llg.E(state) <<std::endl;
    std::cout << "Analytical= " << 0 <<std::endl;
    std::cout << "H_dip_1   = " << -2*param.p/(4*M_PI*pow(dx,3))<<"," <<0<<","<<0<<std::endl;
    std::cout << "H_dip_2   = " <<0<<","<<0<<","<< param.p/(4*M_PI*pow(dx,3)) <<std::endl;
    //-------------------------------------------------------
    m = constant(0.0,mesh.n0,mesh.n1,mesh.n2,3,f64);
    m(0,0,0,0) = 0;
    m(0,0,0,1) = 0;
    m(0,0,0,2) = 1;

    m(1,0,0,0) = 0;
    m(1,0,0,1) = 0;
    m(1,0,0,2) =-1;
    state.m=m;
    vti_writer_atom(state.m, mesh ,(filepath + "/minit2").c_str());
    std::cout << "ENERGY    = " << Llg.E(state) <<std::endl;
    std::cout << "Analytical= " <<  pow(param.p,2)*param.mu0/(4.*M_PI)/pow(dx,3) <<std::endl;
    std::cout << "H_dip_1   = " <<0<<","<<0<<","<<-param.p/(4*M_PI*pow(dx,3)) <<std::endl;
    std::cout << "H_dip_2   = " <<0<<","<<0<<","<< param.p/(4*M_PI*pow(dx,3)) <<std::endl;
    //-------------------------------------------------------

    nx = 1, ny=2 ,nz=1;
    mesh=Mesh(nx,ny,nz,dx,dx,dx);
    m = constant(0.0,mesh.n0,mesh.n1,mesh.n2,3,f64);
    m(0,0,0,0) = 0;
    m(0,0,0,1) = 0;
    m(0,0,0,2) = 1;

    m(0,1,0,0) = 0;
    m(0,1,0,1) = 0;
    m(0,1,0,2) = 1;
    state = State (mesh,param, m);
    vti_writer_atom(state.m, mesh ,(filepath + "/minit").c_str());
  
    llgterm.pop_back();
    llgterm.push_back( llgt_ptr (new ATOMISTIC_DEMAG(mesh)));
    //TODO this leads to compiler error
    //Llg=LLG(state,llgterm);
    LLG Llg2(state,llgterm);
    std::cout << "ENERGY    = " << Llg2.E(state) <<std::endl;
    std::cout << "Analytical= " << - pow(param.p,2)*param.mu0/(4.*M_PI)/pow(dx,3) <<std::endl;//TODO calc on paper, but should be like case 1
    //-------------------------------------------------------
    m = constant(0.0,mesh.n0,mesh.n1,mesh.n2,3,f64);
    m(0,0,0,0) = 1;
    m(0,0,0,1) = 0;
    m(0,0,0,2) = 0;

    m(0,1,0,0) = 0;
    m(0,1,0,1) = 0;
    m(0,1,0,2) = 1;
    state.m=m;
    std::cout << "ENERGY    = " << Llg2.E(state) <<std::endl;
    std::cout << "Analytical= " << 0 <<std::endl;
    return 0;
}