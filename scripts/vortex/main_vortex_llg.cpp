#include "arrayfire.h"
#include "magnum_af.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace af; typedef std::shared_ptr<LLGTerm> llgt_ptr; 

void calc_mean_m(const State& state, const long int n_cells,  std::ostream& myfile){
    myfile << std::setw(12) << state.t << "\t" << afvalue(sum(sum(sum(state.m(span,span,span,0),0),1),2))/n_cells << std::endl;
}

void calc_mean_m(const State& state, const long int n_cells,  std::ostream& myfile, double hzee){
    //TODO test n_cells accuracy
    //myfile << std::setw(12) << state.t << "\t" << afvalue(sum(sum(sum(state.m(span,span,span,0),0),1),2))/n_cells << "\t" << hzee << std::endl;
    array sum_dim3 = sum(sum(sum(state.m,0),1),2);
    myfile << std::setw(12) << state.t << "\t" << afvalue(sum_dim3(span,span,span,0))/n_cells << "\t" << afvalue(sum_dim3(span,span,span,1))/n_cells<< "\t" << afvalue(sum_dim3(span,span,span,2))/n_cells << "\t" << hzee << std::endl;
}

af::array zee_func(State state){
    double field_Tesla = 0;
    double rate = 0.34e6 ; //[T/s]
    double hzee_max = 0.12; //[T]
    if(state.t < hzee_max/rate) field_Tesla = rate *state.t; 
    else if (state.t < 3*hzee_max/rate) field_Tesla = -rate *state.t + 2*hzee_max; 
    else if(state.t < 4*hzee_max/rate) field_Tesla = rate*state.t - 4*hzee_max; 
    else {field_Tesla = 0; std::cout << "WARNING ZEE time out of range" << std::endl;}
    array zee = constant(0.0,state.mesh.n0,state.mesh.n1,state.mesh.n2,3,f64);
    zee(span,span,span,0)=constant(field_Tesla/state.param.mu0 ,state.mesh.n0,state.mesh.n1,state.mesh.n2,1,f64);
    return  zee;
}
  
int main(int argc, char** argv)
{
    std::cout<<"argc"<<argc<<std::endl;
    for (int i=0; i<argc; i++) cout << "Parameter " << i << " was " << argv[i] << "\n";
    std::string filepath(argc>1? argv[1]: "../Data/Testing");
    if(argc>0)filepath.append("/");
    std::cout<<"Writing into path "<<filepath.c_str()<<std::endl;
    setDevice(argc>2? std::stoi(argv[2]):0);
    std::string path_mrelax(argc>3? argv[3]: "");
    info();
    std::cout.precision(24);

    // Parameter initialization
    const int nx = 250, ny=250 ,nz=1;
    //const int nx = 400, ny=400 ,nz=1;
    const double x=1600e-9, y=1600e-9, z=65e-9;//[m] // Physical dimensions
  
    //Generating Objects
    Mesh mesh(nx,ny,nz,x/nx,y/ny,z/nz);
    Param param = Param();
    param.ms    = 1.393e6;//[J/T/m^3] == [Joule/Tesla/meter^3] = 1.75 T/mu_0
    param.A     = 1.5e-11;//[J/m]
    param.alpha = 0.02;
    // Initial magnetic field
    array Ms = constant(0.0,mesh.n0,mesh.n1,mesh.n2,3,f64);
    array m = constant(0.0,mesh.n0,mesh.n1,mesh.n2,3,f64);
    long int n_cells=0;//Number of cells with Ms!=0
    for(int ix=0;ix<mesh.n0;ix++){
        for(int iy=0;iy<mesh.n1;iy++){
            const double rx=double(ix)-mesh.n0/2.;
            const double ry=double(iy)-mesh.n1/2.;
            const double r = sqrt(pow(rx,2)+pow(ry,2));
            if(r<nx/2.){
                for(int iz=0;iz<mesh.n2;iz++){
                    n_cells++;
                }
                Ms(ix,iy,span,span)=param.ms;
                if(r==0.){
                    m(ix,iy,span,2)= 1;
                }
                else{
                    m(ix,iy,span,0)=-ry/r;
                    m(ix,iy,span,1)= rx/r;
                    m(ix,iy,span,2)= sqrt(nx)/r;
                }
            }
        }
    }

    std::cout << "n_cells= " << n_cells << ", should be nx^2*M_PI/4.= " << pow(nx,2)*M_PI/4. << std::endl;

    State state(mesh,param, m);
    state.Ms=Ms;
    vti_writer_micro(state.Ms, mesh ,(filepath + "Ms").c_str());
    calc_mean_m(state,n_cells,std::cout);

    std::vector<llgt_ptr> llgterm;
    llgterm.push_back( llgt_ptr (new DemagSolver(mesh,param)));
    llgterm.push_back( llgt_ptr (new ExchSolver(mesh,param)));
    LLG Llg(state,llgterm);

    //obtaining relaxed magnetization
    if(!exists (path_mrelax)){
        std::cout << "mrelax.vti not found, starting relaxation" << std::endl;
        vti_writer_micro(state.m, mesh ,(filepath + "minit_nonnormalized").c_str());
        state.m=renormalize_handle_zero_values(state.m);
        vti_writer_micro(state.m, mesh ,(filepath + "minit_renorm").c_str());
        timer t = af::timer::start();
        double E_prev=1e20;
        double E=0;
        while (fabs((E_prev-E)/E_prev) > 1e-7){
            state.m=Llg.llgstep(state);
            if( state.steps % 1000 == 0){
                E_prev=E;
                E=Llg.E(state);
                std::cout << "step " << state.steps << " time= " << state.t << " E= " << E << " fabs= " << fabs((E_prev-E)/E_prev) << std::endl;
                vti_writer_micro(state.m, mesh ,(filepath + "m_current_relaxation"+std::to_string(state.steps)).c_str());
            }
        }
        std::cout<<"timerelax [af-s]: "<< af::timer::stop(t) <<std::endl;
        vti_writer_micro(state.m, mesh ,(filepath + "mrelax").c_str());
        std::cout << "state.t= "<< state.t << std::endl;
        state.t=0;//setting new zero time
    }
    else{
        std::cout << "found mrelax. loading magnetization" << std::endl;
        vti_reader(m, mesh, filepath+"mrelax.vti");
        state.m=m;
    }
    vti_writer_micro(state.m, mesh ,(filepath + "todel_check_mrelax").c_str());

    std::ofstream stream;
    stream.precision(12);
    stream.open ((filepath + "m.dat").c_str());
    stream << "# t	<mx>" << std::endl;
    calc_mean_m(state,n_cells,stream);

    timer t_hys = af::timer::start();
    double rate = 0.34e6 ; //[T/s]
    double hzee_max = 0.12; //[T]
    Llg.Fieldterms.push_back( llgt_ptr (new Zee(&zee_func))); //Rate in T/s
    while (state.t < 4* hzee_max/rate){
         state.m=Llg.llgstep(state);
         calc_mean_m(state,n_cells,stream,afvalue(Llg.Fieldterms[2]->h(state)(0,0,0,0)));
         if( state.steps % 2000 == 0){
             vti_writer_micro(state.m, mesh ,(filepath + "m_hysteresis_"+std::to_string(state.steps)).c_str());
         }
         //std::cout << state.t << "\t" << afvalue(sum(sum(sum(sum(Llg.Fieldterms[2]->h(state),3),2),1),0))<< std::endl;
    }
    stream.close();
    std::cout<<"time full hysteresis [af-s]: "<< af::timer::stop(t_hys) <<std::endl;
    return 0;
}
