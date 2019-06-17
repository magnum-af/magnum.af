#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace af;
typedef std::shared_ptr<LLGTerm> llgt_ptr;

void calc_mean_m(const State& state, std::ostream& myfile, double hzee){
    const array sum_dim3 = sum(sum(sum(state.m,0),1),2);
    const int ncells = state.mesh.n0 * state.mesh.n1 * state.mesh.n2;
    myfile << std::setw(12) << state.t << "\t" << afvalue(sum_dim3(span,span,span,0))/ncells << "\t" << afvalue(sum_dim3(span,span,span,1))/ncells<< "\t" << afvalue(sum_dim3(span,span,span,2))/ncells << "\t" << hzee << std::endl;
}

const double hzee_max = 2; //[T]
const double simtime = 100e-9;
const double rate = hzee_max/simtime; //[T/s]

af::array zee_func(State state){
    double field_Tesla = 0;
    field_Tesla = rate *state.t;
    array zee = constant(0.0,state.mesh.n0,state.mesh.n1,state.mesh.n2,3,f64);
    zee(span,span,span,0)=constant(field_Tesla/state.constants::mu0 ,state.mesh.n0,state.mesh.n1,state.mesh.n2,1,f64);
    return  zee;
}

int main(int argc, char** argv)
{

    std::cout<<"argc = "<<argc<<std::endl;
     for (int i=0; i<argc; i++)
          cout << "Parameter " << i << " was " << argv[i] << "\n";

    std::string filepath(argc>1? argv[1]: "../Data/skyrmion_stoch");
    if(argc>0)filepath.append("/");
    std::cout<<"Writing into path "<<filepath.c_str()<<std::endl;

    setDevice(argc>2? std::stoi(argv[2]):0);
    info();

    // Parameter initialization
    const int nx = 90, nz=1;
    const double dx=1.0e-9;
    const double dz=0.6e-9;

    //Generating Objects
    Mesh mesh(nx,nx,nz,dx,dx,dz);
    Material material = Material();
    state.Ms    = 580000;
    material.A     = 15e-12;
    material.alpha = 1;
    material.D=3e-3;
    material.Ku1=0.6e6;

     // Initial magnetic field
     array m = constant(0.0,mesh.n0,mesh.n1,mesh.n2,3,f64);
     m(span,span,span,2) = -1;
     for(int ix=0;ix<mesh.n0;ix++){
         for(int iy=0;iy<mesh.n1;iy++){
             const double rx=double(ix)-mesh.n0/2.;
             const double ry=double(iy)-mesh.n1/2.;
             const double r = sqrt(pow(rx,2)+pow(ry,2));
             if(r>nx/4.) m(ix,iy,span,2)=1.;
         }
     }

    State state(mesh,material, m);
    vti_writer_atom(state.m, mesh ,(filepath + "minit").c_str());

    // Relax
    std::vector<llgt_ptr> llgterm;
    llgterm.push_back( llgt_ptr (new DemagField(mesh,material)));
    llgterm.push_back( llgt_ptr (new ExchangeField(mesh,material)));
    llgterm.push_back( llgt_ptr (new DmiField(mesh,material)));
    llgterm.push_back( llgt_ptr (new UniaxialAnisotropyField(mesh,material)));

    LLG Llg(state,llgterm);

    timer t = af::timer::start();
    double E_prev=1e20;
    while (fabs((E_prev-Llg.E(state))/E_prev) > 1e-8){
        E_prev=Llg.E(state);
        for ( int i = 0; i<100; i++){
            state.m=Llg.step(state);
        }
        if( state.steps % 1000 == 0) std::cout << "step " << state.steps << " rdiff= " << fabs((E_prev-Llg.E(state))/E_prev) << std::endl;
    }
    std::cout << "time =" << state.t << " [s], E = " << Llg.E(state) << "[J]" << std::endl;
    std::cout<<"timerelax [af-s]: "<< af::timer::stop(t) << ", steps = " << state.steps << std::endl;
    vti_writer_micro(state.m, mesh ,(filepath + "relax").c_str());

    // Hysteresis
    std::ofstream stream;
    stream.precision(12);
    stream.open ((filepath + "m.dat").c_str());
    stream << "# t	<mx>    <my>    <mz>    hzee" << std::endl;

    timer t_hys = af::timer::start();
    Llg.Fieldterms.push_back( llgt_ptr (new ExternalField(&zee_func))); //Rate in T/s
    while (state.t <  simtime){
         state.m=Llg.step(state);
         if( state.steps % 1000 == 0){
             calc_mean_m(state, stream, afvalue(Llg.Fieldterms[Llg.Fieldterms.size()-1]->h(state)(0,0,0,0)));
             vti_writer_micro(state.m, mesh ,(filepath + "m_hysteresis_"+std::to_string(state.steps)).c_str());
         }
    }
    stream.close();
    std::cout<<"time full hysteresis [af-s]: "<< af::timer::stop(t_hys) <<std::endl;
    return 0;
}
