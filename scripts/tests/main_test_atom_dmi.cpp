#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace magnumafcpp;

using namespace af; typedef std::shared_ptr<LLGTerm> llgt_ptr;

bool compare(double a, double b){
    //std::cout << "COM:"<< a <<", " << b <<", "<<fabs(a-b)/fabs(a+b)<<std::endl;
    if(a == 0 && b == 0) return false;
    if(fabs(a-b)/fabs(a+b)<1e-30) return false; //TODO fore some values of D_atom and p, this is too low, use e.g. 1e-15
    else return true;
}

int main(int argc, char** argv)
{
    info();
    std::cout.precision(32);
    std::string filepath(argc>0? argv[1]: "../Data/Testing/");
    int nx = 2, ny=1 , nz=1;//nz=5 -> lz=(5-1)*dx
    //const double dx=1;
    const double dx=2.715e-10;

    //Generating Objects
    Mesh mesh(nx, ny, nz, dx, dx, dx);
    Material material = Material();
    material.D_atom=1e23;
    material.p=1e-22;
    //material.p=9.274009994e-24;

    //-------------------------------------------------------
    array m = constant(0.0, mesh.n0, mesh.n1, mesh.n2, 3, f64);
    m(0, 0, 0, 0) = 0;
    m(0, 0, 0, 1) = 0;
    m(0, 0, 0, 2) = 1;

    m(1, 0, 0, 0) = 0;
    m(1, 0, 0, 1) = 0;
    m(1, 0, 0, 2) = 1;
    State state(mesh, material, m);
    vti_writer_atom(state.m, mesh , (filepath + "/minit").c_str());

    std::vector<llgt_ptr> llgterm;
    llgterm.push_back( llgt_ptr (new AtomisticDmiField(mesh, material)));
    LLG Llg(state, llgterm);
    if(compare(Llg.E(state), 0)) std::cout <<"!!! TEST 1 FAILED !!!"<< std::endl;
    if(compare(afvalue(Llg.Fieldterms[0]->h(state)(0, 0, 0, 0)), -material.D_atom/constants::mu0/material.p)) std::cout <<"!!! TEST 1 H1x FAILED !!!"<< std::endl;
    if(compare(afvalue(Llg.Fieldterms[0]->h(state)(0, 0, 0, 1)), 0)) std::cout <<"!!! TEST 1 H1y FAILED !!!"<< std::endl;
    if(compare(afvalue(Llg.Fieldterms[0]->h(state)(0, 0, 0, 2)), 0)) std::cout <<"!!! TEST 1 H1z FAILED !!!"<< std::endl;
    if(compare(afvalue(Llg.Fieldterms[0]->h(state)(1, 0, 0, 0)), material.D_atom/constants::mu0/material.p)) std::cout <<"!!! TEST 1 H2x FAILED !!!"<< std::endl;
    if(compare(afvalue(Llg.Fieldterms[0]->h(state)(1, 0, 0, 1)), 0)) std::cout <<"!!! TEST 1 H2y FAILED !!!"<< std::endl;
    if(compare(afvalue(Llg.Fieldterms[0]->h(state)(1, 0, 0, 2)), 0)) std::cout <<"!!! TEST 1 H2z FAILED !!!"<< std::endl;
    //std::cout << "ENERGY    = " << Llg.E(state) <<std::endl;
    //std::cout << "Analytical= " << analytical <<std::endl;
    //std::cout << "H_exch_1   = " <<0<<", "<<0<<", "<< material.J_atom/constants::mu0/material.p <<std::endl;
    //std::cout << "H_exch_2   = " <<0<<", "<<0<<", "<< material.J_atom/constants::mu0/material.p <<std::endl;

//    //-------------------------------------------------------
    m = constant(0.0, mesh.n0, mesh.n1, mesh.n2, 3, f64);
    m(0, 0, 0, 0) = 0;
    m(0, 0, 0, 1) = 0;
    m(0, 0, 0, 2) = 1;

    m(1, 0, 0, 0) = 1;
    m(1, 0, 0, 1) = 0;
    m(1, 0, 0, 2) = 0;

    state.m=m;
    //vti_writer_atom(state.m, mesh , (filepath + "/minit2").c_str());
    if(compare(Llg.E(state), -material.D_atom)) std::cout <<"!!! TEST 2 FAILED !!!"<< std::endl;
    if(compare(afvalue(Llg.Fieldterms[0]->h(state)(0, 0, 0, 0)), 0)) std::cout <<"!!! TEST 2 H1x FAILED !!!"<< std::endl;
    if(compare(afvalue(Llg.Fieldterms[0]->h(state)(0, 0, 0, 1)), 0)) std::cout <<"!!! TEST 2 H1y FAILED !!!"<< std::endl;
    if(compare(afvalue(Llg.Fieldterms[0]->h(state)(0, 0, 0, 2)), material.D_atom/constants::mu0/material.p)) std::cout <<"!!! TEST 2 H1z FAILED !!!"<< std::endl;
    if(compare(afvalue(Llg.Fieldterms[0]->h(state)(1, 0, 0, 0)), material.D_atom/constants::mu0/material.p)) std::cout <<"!!! TEST 2 FAILED !!!"<< std::endl;
    if(compare(afvalue(Llg.Fieldterms[0]->h(state)(1, 0, 0, 1)), 0)) std::cout <<"!!! TEST 2 FAILED !!!"<< std::endl;
    if(compare(afvalue(Llg.Fieldterms[0]->h(state)(1, 0, 0, 2)), 0)) std::cout <<"!!! TEST 2 FAILED !!!"<< std::endl;
    //std::cout << "ENERGY    = " << Llg.E(state) <<std::endl;
    //std::cout << "Analytical= " << 0 <<std::endl;
    //std::cout << "H_exch_1   = " <<material.J_atom/constants::mu0/material.p<<", "<<0<<", "<<0  <<std::endl;
    //std::cout << "H_exch_2   = " <<0<<", "<<0<<", "<< material.J_atom/constants::mu0/material.p <<std::endl;
    //-------------------------------------------------------
////    m = constant(0.0, mesh.n0, mesh.n1, mesh.n2, 3, f64);
////    m(0, 0, 0, 0) = 0;
////    m(0, 0, 0, 1) = 0;
////    m(0, 0, 0, 2) = 1;
////
////    m(1, 0, 0, 0) = 0;
////    m(1, 0, 0, 1) = 0;
////    m(1, 0, 0, 2) =-1;
////    state.m=m;
////    if(compare(Llg.E(state), material.J_atom)) std::cout <<"!!! TEST 3 FAILED !!!"<< std::endl;
////    if(compare(afvalue(Llg.Fieldterms[0]->h(state)(0, 0, 0, 0)), 0)) std::cout <<"!!! TEST 3 FAILED !!!"<< std::endl;
////    if(compare(afvalue(Llg.Fieldterms[0]->h(state)(0, 0, 0, 1)), 0)) std::cout <<"!!! TEST 3 FAILED !!!"<< std::endl;
////    if(compare(afvalue(Llg.Fieldterms[0]->h(state)(0, 0, 0, 2)), -material.J_atom/constants::mu0/material.p)) std::cout <<"!!! TEST 3 FAILED !!!"<< std::endl;
////    if(compare(afvalue(Llg.Fieldterms[0]->h(state)(1, 0, 0, 0)), 0)) std::cout <<"!!! TEST 3 FAILED !!!"<< std::endl;
////    if(compare(afvalue(Llg.Fieldterms[0]->h(state)(1, 0, 0, 1)), 0)) std::cout <<"!!! TEST 3 FAILED !!!"<< std::endl;
////    if(compare(afvalue(Llg.Fieldterms[0]->h(state)(1, 0, 0, 2)), material.J_atom/constants::mu0/material.p)) std::cout <<"!!! TEST 3 FAILED !!!"<< std::endl;
////    //std::cout << "ENERGY    = " << Llg.E(state) <<std::endl;
////    //std::cout << "Analytical= " << material.J_atom <<std::endl;
////    //std::cout << "H_exch_1   = " <<0<<", "<<0<<", "<<-material.J_atom/constants::mu0/material.p <<std::endl;
////    //std::cout << "H_exch_2   = " <<0<<", "<<0<<", "<< material.J_atom/constants::mu0/material.p <<std::endl;
////    //std::cout << "Analytical= " << afvalue(Llg.Fieldterms[0]->h(state)(0, 0, 0, 2))<<std::endl;
////    //af::print("LLG", Llg.Fieldterms[0]->h(state));
//////    //-------------------------------------------------------
//////
    nx = 1, ny=1 , nz=2;
    mesh=Mesh(nx, ny, nz, dx, dx, dx);
    m = constant(0.0, mesh.n0, mesh.n1, mesh.n2, 3, f64);
    m(0, 0, 0, 0) = 0;
    m(0, 0, 0, 1) = 0;
    m(0, 0, 0, 2) = 1;

    m(0, 0, 1, 0) = 1;
    m(0, 0, 1, 1) = 0;
    m(0, 0, 1, 2) = 0;
    state = State (mesh, material, m);

    llgterm.pop_back();
    llgterm.push_back( llgt_ptr (new AtomisticDmiField(mesh, material)));
    //TODO this leads to compiler error
    //Llg=LLG(state, llgterm);
    LLG Llg2(state, llgterm);
    if(compare(Llg2.E(state), 0)) std::cout <<"!!! TEST 4 FAILED !!!"<< std::endl;
    if(compare(afvalue(Llg2.Fieldterms[0]->h(state)(0, 0, 0, 0)), 0)) std::cout <<"!!! TEST 4 FAILED !!!"<< std::endl;
    if(compare(afvalue(Llg2.Fieldterms[0]->h(state)(0, 0, 0, 1)), 0)) std::cout <<"!!! TEST 4 FAILED !!!"<< std::endl;
    if(compare(afvalue(Llg2.Fieldterms[0]->h(state)(0, 0, 0, 2)), 0)) std::cout <<"!!! TEST 4 FAILED !!!"<< std::endl;
    if(compare(afvalue(Llg2.Fieldterms[0]->h(state)(0, 0, 1, 0)), 0)) std::cout <<"!!! TEST 4 FAILED !!!"<< std::endl;
    if(compare(afvalue(Llg2.Fieldterms[0]->h(state)(0, 0, 1, 1)), 0)) std::cout <<"!!! TEST 4 FAILED !!!"<< std::endl;
    if(compare(afvalue(Llg2.Fieldterms[0]->h(state)(0, 0, 1, 2)), 0)) std::cout <<"!!! TEST 4 FAILED !!!"<< std::endl;
//////    //-------------------------------------------------------
//////    m = constant(0.0, mesh.n0, mesh.n1, mesh.n2, 3, f64);
//////    m(0, 0, 0, 0) = 1;
//////    m(0, 0, 0, 1) = 0;
//////    m(0, 0, 0, 2) = 0;
//////
//////    m(0, 1, 0, 0) = 0;
//////    m(0, 1, 0, 1) = 0;
//////    m(0, 1, 0, 2) = 1;
//////    state.m=m;
//////    std::cout << "ENERGY    = " << Llg2.E(state) <<std::endl;
//////    std::cout << "Analytical= " << 0 <<std::endl;
    return 0;
}
