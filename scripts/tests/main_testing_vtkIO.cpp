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
using namespace af; typedef std::shared_ptr<LLGTerm> llgt_ptr; 
void calcm(State state, std::ostream& myfile);
int main(int argc, char** argv)
{
  std::cout<<"argc"<<argc<<std::endl;
   for (int i=0; i<argc; i++)
        cout << "Parameter " << i << " was " << argv[i] << "\n";
  
  std::string filepath(argc>1? argv[1]: "../Data/Testing");
  if(argc>0)filepath.append("/");
  std::cout<<"Writing into path "<<filepath.c_str()<<std::endl;

  setDevice(argc>2? std::stoi(argv[2]):0);
  //if(argc>1) setDevice(std::stoi(argv[2]));
  info();

  // Parameter initialization
  const double x=5.e-7, y=1.25e-7, z=3.e-9;
  const int nx = 5, ny=4 ,nz=2;


  //Simulation Parameters
  //double hmax = 3.5e-10;
  //double hmin = 1.0e-15;
  //double atol = 1e-8;
  //double rtol = atol;
  
  //Generating Objects
  Mesh mesh(nx,ny,nz,x/nx,y/ny,z/nz);
  Material material = Material();
  material.ms    = 8e5;
  material.A     = 1.3e-11;
  material.alpha = 1;
  material.afsync  = false;

  // Initial magnetic field
  array m = constant(0.0,mesh.n0,mesh.n1,mesh.n2,3,f64);
  m(seq(1,end-1),span,span,0) = constant(1.0,mesh.n0-2,mesh.n1,mesh.n2,1,f64);
  m(0,span,span,1 ) = constant(1.0,1,mesh.n1,mesh.n2,1,f64);
  m(-1,span,span,1) = constant(1.0,1,mesh.n1,mesh.n2,1,f64);
  m=iota(dim4(nx,ny,nz,3),dim4(1,1,1,1),f64);
  print("A",m);
  State state(mesh,material, m);
  vti_writer_micro(state.m, mesh ,(filepath + "minit").c_str());
  array A = array();
  Mesh newmesh = Mesh(0,0,0,0,0,0);
  vti_reader(A,newmesh,"/home/pth/git/magnum.af/Data/Testing/minit.vti");

  vti_writer_atom(A, newmesh ,(filepath + "aminit").c_str());
  print("A",A);
  std::cout<<newmesh.n0<<"  "<<newmesh.n1<<"  "<<newmesh.n2<<"  "<<newmesh.dx<<"  "<<newmesh.dy<<"  "<<newmesh.dz<<"  "<<std::endl;
  array B = array();
  Mesh bmesh = Mesh(0,0,0,0,0,0);

  vti_reader(B,bmesh,"/home/pth/git/magnum.af/Data/Testing/aminit.vti");
 
  print("B",B);
  std::cout<<bmesh.n0<<"  "<<bmesh.n1<<"  "<<bmesh.n2<<"  "<<bmesh.dx<<"  "<<bmesh.dy<<"  "<<bmesh.dz<<"  "<<std::endl;

  vtr_writer(B, bmesh ,(filepath + "minit").c_str());
  array C = array();
  Mesh cmesh = Mesh(0,0,0,0,0,0);
  vtr_reader(C,cmesh,"/home/pth/git/magnum.af/Data/Testing/minit.vtr");
  print("C",C);
  std::cout<<cmesh.n0<<"  "<<cmesh.n1<<"  "<<cmesh.n2<<"  "<<cmesh.dx<<"  "<<cmesh.dy<<"  "<<cmesh.dz<<"  "<<std::endl;
//  std::vector<llgt_ptr> llgterm;
//  llgterm.push_back( llgt_ptr (new DemagField(mesh,material)));
//  llgterm.push_back( llgt_ptr (new ExchangeField(mesh,material)));
//  LLG Llg(state,llgterm);
//  //LLG Llg(state,atol,rtol,hmax,hmin,llgterm);
//
//  std::ofstream stream;
//  stream.precision(12);
//  stream.open ((filepath + "m.dat").c_str());
//  
//  timer t = af::timer::start();
//  while (state.t < 1.e-9){
//    state.m=Llg.llgstep(state);
//    calcm(state,stream);
//  }
//  std::cout<<"Energy of relaxed state = "<<Llg.E(state)<<"\n"<<std::endl;
//  double timerelax= af::timer::stop(t);
//  af_to_vti(state.m, mesh ,(filepath + "relax").c_str());
//
//  std::cout<<"timerelax [af-s]: "<< timerelax << " for "<<Llg.counter_accepted+Llg.counter_reject<<" steps, thereof "<< Llg.counter_accepted << " Steps accepted, "<< Llg.counter_reject<< " Steps rejected" << std::endl;
//
//  // Prepare switch
//  array zeeswitch = constant(0.0,1,1,1,3,f64);
//  zeeswitch(0,0,0,0)=-24.6e-3/constants::mu0;
//  zeeswitch(0,0,0,1)=+4.3e-3/constants::mu0;
//  zeeswitch(0,0,0,2)=0.0;
//  zeeswitch = tile(zeeswitch,mesh.n0,mesh.n1,mesh.n2);
//  llgterm.push_back( llgt_ptr (new Zee(zeeswitch,mesh,material)));
//  Llg.Fieldterms=llgterm;
//  //TODO remove state0 in LLG!
//  Llg.state0.material.alpha=0.02;
//
//  while (state.t < 2.e-9){
//    state.m=Llg.llgstep(state);
//    calcm(state,stream);
//  }
//  af_to_vti(state.m, mesh ,(filepath + "2ns").c_str());
//  stream.close();
//  Llg.print_cpu_time(std::cout); 
  return 0;
}
void calcm(State state, std::ostream& myfile){
  myfile << std::setw(12) << state.t << "\t" <<meani(state.m,0)<< "\t" <<meani(state.m,1)<< "\t" <<meani(state.m,2)<< "\t" << std::endl;
}
