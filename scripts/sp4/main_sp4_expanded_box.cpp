#include "arrayfire.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <string>
#include <memory>
#include "llg.hpp"
#include "micro_exch.hpp"
#include "micro_demag.hpp"
#include "zee.hpp"
#include "atomistic_demag.hpp"
#include "atomistic_exchange.hpp"
#include "atomistic_anisotropy.hpp"
#include "atomistic_dmi.hpp"
#include "vtk_IO.hpp"
#include "string.hpp"
#include "stochastic_llg.hpp"

using namespace magnumaf;

using namespace af; typedef std::shared_ptr<LLGTerm> llgt_ptr;
//void calcm(State state, std::ostream& myfile);
void calcm(State state, std::ostream& myfile, int nx, int ny, int nz, int spnx, int spny, int spnz);
void calcm_autodetect(State state, std::ostream& myfile){
    //TODOsum(state.m(span, span, span, 0), 0)
    myfile << std::setw(12) << state.t << "\t" <<meani(state.m, 0)<< "\t" <<meani(state.m, 1)<< "\t" <<meani(state.m, 2)<< "\t" << std::endl;
}

int main(int argc, char** argv)
{
    std::cout<<"argc"<<argc<<std::endl;
    for (int i=0; i<argc; i++) std::cout << "Parameter " << i << " was " << argv[i] << "\n";
    std::string filepath(argc>1? argv[1]: "../Data/Testing");
    if(argc>0)filepath.append("/");
    std::cout<<"Writing into path "<<filepath.c_str()<<std::endl;
    setDevice(argc>2? std::stoi(argv[2]):0);
    info();

//    //TODO  Include in python testing
//    //TEST RENORMALIZE
//    array test_renorm=constant(0.0, 6, 1, 1, 3, f64);
//    //nx=0: (1, 0, 0)
//    test_renorm(0, 0, 0, 0)=1;
//    test_renorm(0, 0, 0, 1)=0;
//    test_renorm(0, 0, 0, 2)=0;
//    //nx=1: (1, 0, 0)
//    test_renorm(1, 0, 0, 0)=0;
//    test_renorm(1, 0, 0, 1)=0;
//    test_renorm(1, 0, 0, 2)=0;
//    //nx=2: (1, 0, 0)
//    test_renorm(2, 0, 0, 0)=0;
//    test_renorm(2, 0, 0, 1)=1;
//    test_renorm(2, 0, 0, 2)=0;
//    //nx=3: (1, 0, 0)
//    test_renorm(3, 0, 0, 0)=0;
//    test_renorm(3, 0, 0, 1)=0;
//    test_renorm(3, 0, 0, 2)=1;
//    //nx=4: (1, 0, 0)
//    test_renorm(4, 0, 0, 0)=1;
//    test_renorm(4, 0, 0, 1)=1;
//    test_renorm(4, 0, 0, 2)=1;
//    //nx=5: (1, 0, 0)
//    test_renorm(5, 0, 0, 0)=2;
//    test_renorm(5, 0, 0, 1)=2;
//    test_renorm(5, 0, 0, 2)=2;
//    print ("test_renorm", test_renorm);
//    print ("renormalized test_renorm", renormalize(test_renorm));
//    print ("renormalized handle zero  test_renorm", renormalize_handle_zero_values(test_renorm));
//    //END TEST RENORMALIZE


    // Parameter initialization
    //const double x=5.e-7, y=1.25e-7, z=3.e-9;
    //const int nx = 100, ny=25 , nz=1; //NOTE This with CASE 1 yields same results as backups
    const int nx = 120, ny=45 , nz=3; //Total Box dimensions
    const int spnx = 100, spny=25 , spnz=1;// Box dimensions without vacuum

    //Generating Objects
    Mesh mesh(nx, ny, nz, 5.e-7/100, 1.25e-7/25, 3.e-9);
    Material material = Material();
    state.Ms    = 8e5;
    material.A     = 1.3e-11;
    material.alpha = 1;

    // Initial magnetic field
    //CASE 1
    //array m = constant(0.0, mesh.n0, mesh.n1, mesh.n2, 3, f64);
    //m(seq(1, end-1), span, span, 0) = constant(1.0, mesh.n0-2, mesh.n1, mesh.n2, 1, f64);
    //m(0, span, span, 1 ) = constant(1.0, 1, mesh.n1, mesh.n2, 1, f64);
    //m(-1, span, span, 1) = constant(1.0, 1, mesh.n1, mesh.n2, 1, f64);
    //CASE 1

    //CASE 2
    array m = constant(0.0, mesh.n0, mesh.n1, mesh.n2, 3, f64);
    m(seq((nx-spnx)/2+1, end-(nx-spnx)/2-1), seq((ny-spny)/2, end-(ny-spny)/2), seq((nz-spnz)/2, end-(nz-spnz)/2), 0) = constant(1.0, spnx-2, spny, spnz, 1, f64);
    m((ny-spny)/2, seq((ny-spny)/2, end-(ny-spny)/2), seq((nz-spnz)/2, end-(nz-spnz)/2), 1 ) = constant(1.0, 1, spny, spnz, 1, f64);
    m(-(ny-spny)/2-1, seq((ny-spny)/2, end-(ny-spny)/2), seq((nz-spnz)/2, end-(nz-spnz)/2), 1) = constant(1.0, 1, spny, spnz, 1, f64);
    //std::cout << m(seq((nx-spnx)/2+1, end-(nx-spnx)/2-1), seq((ny-spny)/2, end-(ny-spny)/2), seq((nz-spnz)/2, end-(nz-spnz)/2), 0).dims() << std::endl;
    //std::cout << constant(1.0, mesh.n0-2, mesh.n1, mesh.n2, 1, f64).dims() << std::endl;
    //std::cout << m(-(ny-spny)/2-1, seq((ny-spny)/2, end-(ny-spny)/2), seq((nz-spnz)/2, end-(nz-spnz)/2), 1).dims() << std::endl;
    //print ("case 2 m:", m);
    //CASE 2

    State state(mesh, material, m);
    vti_writer_micro(state.m, mesh , (filepath + "minit").c_str());

    //CASE 1
    //state.Ms=constant(state.Ms, state.mesh.dims, f64);
    //CASE 1

    //CASE 2
    state.Ms =constant(0.0, state.mesh.dims, f64);
    state.Ms(seq((nx-spnx)/2, end-(nx-spnx)/2), seq((ny-spny)/2, end-(ny-spny)/2), seq((nz-spnz)/2, end-(nz-spnz)/2), span) = constant(state.Ms, spnx, spny, spnz, 3, f64);
    //CASE 2

    vti_writer_micro(state.Ms, mesh , (filepath + "Ms").c_str());

    //testing MS
    //std::cout << "is_empty: "<< state.Ms.isempty()<<std::endl;
    //state.Ms=constant(state.Ms, state.mesh.dims, f64);
    //std::cout << "is_empty: "<< state.Ms.isempty()<<std::endl;
    //(state.Ms.isempty()? std::cout << "TRUE" << true <<std::endl : std::cout << "flase" << false <<std::endl);
    //if (state.Ms.isempty()){std::cout << "state.MS.isempts()"<<std::endl;}
    //else {std::cout << "state.MS.is NOT empty()"<<std::endl;}

    //mesh=Mesh(4, 4, 4, x/nx, y/ny, z/nz);
    //m=constant(0, mesh.dims, f64);
    //m(span, span, span, 0) = constant(1.0, mesh.n0 , mesh.n1, mesh.n2, 1, f64);
    //print("m", m);

    //array Ms=constant(1., mesh.dims, f64);
    //Ms(0, 0, 0, span)=constant(0. , 1, 1, 1, 3, f64);
    //print("Ms", Ms);

    //array Div = m/Ms;
    //print("Div", Div);
    ////print("isNaNDiv", isNaN(Div));
    //print("isinfDiv", isInf(Div));
    //print("!isinfDiv", !isInf(Div));
    //print("Div", Div);
    ////replace(Div, !isInf(Div), 0);
    ////replace(Div, !isNaN(Div), 0);
    //replace(Div, Ms!=0, 0);
    //print("Div", Div);

    std::vector<llgt_ptr> llgterm;
    llgterm.push_back( llgt_ptr (new DemagField(mesh, material)));
    llgterm.push_back( llgt_ptr (new ExchangeField(mesh, material)));
    LLG Llg(state, llgterm);

    std::ofstream stream;
    stream.precision(12);
    stream.open ((filepath + "m.dat").c_str());

    timer t = af::timer::start();
    //for (int i = 0; i<5; i++){
    while (state.t < 5.e-10){
      state.m=Llg.step(state);
      calcm(state, std::cout, nx, ny, nz, spnx, spny, spnz);
      calcm(state, stream   , nx, ny, nz, spnx, spny, spnz);
    }
    double timerelax= af::timer::stop(t);
    std::cout<<"timerelax [af-s]: "<< timerelax <<std::endl;
    vti_writer_micro(state.m, mesh , (filepath + "relax").c_str());

    // Prepare switch
    array zeeswitch = constant(0.0, 1, 1, 1, 3, f64);
    zeeswitch(0, 0, 0, 0)=-24.6e-3/constants::mu0;
    zeeswitch(0, 0, 0, 1)=+4.3e-3/constants::mu0;
    zeeswitch(0, 0, 0, 2)=0.0;
    zeeswitch = tile(zeeswitch, mesh.n0, mesh.n1, mesh.n2);
    llgterm.push_back( llgt_ptr (new ExternalField(zeeswitch, mesh, material)));
    Llg.Fieldterms=llgterm;
    Llg.state0.material.alpha=0.02;

    while (state.t < 1.5e-9){
      state.m=Llg.step(state);
      calcm(state, std::cout, nx, ny, nz, spnx, spny, spnz);
      calcm(state, stream   , nx, ny, nz, spnx, spny, spnz);
    }
    vti_writer_micro(state.m, mesh , (filepath + "2ns").c_str());
    stream.close();
    return 0;
}
void calcm(State state, std::ostream& myfile, int nx, int ny, int nz, int spnx, int spny, int spnz){
    myfile << std::setw(12) << state.t << "\t" <<meani(state.m(seq((nx-spnx)/2, end-(nx-spnx)/2), seq((ny-spny)/2, end-(ny-spny)/2), seq((nz-spnz)/2, end-(nz-spnz)/2)), 0)<< "\t" <<meani(state.m(seq((nx-spnx)/2, end-(nx-spnx)/2), seq((ny-spny)/2, end-(ny-spny)/2), seq((nz-spnz)/2, end-(nz-spnz)/2)), 1)<< "\t" <<meani(state.m(seq((nx-spnx)/2, end-(nx-spnx)/2), seq((ny-spny)/2, end-(ny-spny)/2), seq((nz-spnz)/2, end-(nz-spnz)/2)), 2)<< "\t" << std::endl;
}
