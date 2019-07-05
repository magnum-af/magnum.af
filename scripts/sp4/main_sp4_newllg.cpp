#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace magnumaf;


using namespace af; typedef std::shared_ptr<LLGTerm> llgt_ptr;

void calcm(State state, std::ostream& myfile);

int main(int argc, char** argv)
{
     std::cout<<"argc"<<argc<<std::endl;
     for (int i=0; i<argc; i++)
          cout << "Parameter " << i << " was " << argv[i] << "\n";

    std::string filepath(argc>1? argv[1]: "../../Data/Testing");
    if(argc>0)filepath.append("/");
    std::cout<<"Writing into path "<<filepath.c_str()<<std::endl;

    setDevice(argc>2? std::stoi(argv[2]):0);
    info();

    // Parameter initialization
    const double x=5.e-7, y=1.25e-7, z=3.e-9;
    const int nx = 100, ny=25 , nz=1;

    //Generating Objects
    Mesh mesh(nx, ny, nz, x/nx, y/ny, z/nz);
    Material material = Material();
    state.Ms    = 8e5;
    material.A     = 1.3e-11;
    material.alpha = 1;
    state.material.afsync  = false;
    material.T  = 300;

    // Initial magnetic field
    array m = constant(0.0, mesh.n0, mesh.n1, mesh.n2, 3, f64);
    m(seq(1, end-1), span, span, 0) = constant(1.0, mesh.n0-2, mesh.n1, mesh.n2, 1, f64);
    m(0, span, span, 1 ) = constant(1.0, 1, mesh.n1, mesh.n2, 1, f64);
    m(-1, span, span, 1) = constant(1.0, 1, mesh.n1, mesh.n2, 1, f64);
    State state(mesh, material, m);
    vti_writer_micro(state.m, mesh , (filepath + "minit").c_str());

    LLGIntegrator Llg = LLGIntegrator("RKF45");
    Llg.llgterms.push_back( LlgTerm (new DemagField(mesh, material)));
    Llg.llgterms.push_back( LlgTerm (new ExchangeField(mesh, material)));

    std::ofstream stream;
    stream.precision(12);
    stream.open ((filepath + "m.dat").c_str());

    // Relax
    timer t = af::timer::start();
    double time=0;
    while (state.t < 5.e-10){
        timer t2 = af::timer::start();
        Llg.step(state);
        time+=af::timer::stop(t2);
        calcm(state, stream);
    }
    std::cout<<"timerelax              [af-s]: "<< af::timer::stop(t) <<std::endl;
    std::cout<<"time                   [af-s]: "<< time <<std::endl;
    std::cout<<"Llg.get_time_allsteps  [af-s]: "<< Llg.get_time_allsteps() <<std::endl;
    std::cout<<"Llg.get_time_heff      [af-s]: "<< Llg.get_time_heff() <<std::endl;
    vti_writer_micro(state.m, mesh , (filepath + "relax").c_str());

    // Prepare switch
    array zeeswitch = constant(0.0, 1, 1, 1, 3, f64);
    zeeswitch(0, 0, 0, 0)=-24.6e-3/constants::mu0;
    zeeswitch(0, 0, 0, 1)=+4.3e-3/constants::mu0;
    zeeswitch(0, 0, 0, 2)=0.0;
    zeeswitch = tile(zeeswitch, mesh.n0, mesh.n1, mesh.n2);
    Llg.llgterms.push_back( LlgTerm (new ExternalField(zeeswitch)));
    state.material.alpha=0.02;

    // Switch
    t = af::timer::start();
    time=0;
    while (state.t < 1.5e-9){
        timer t2 = af::timer::start();
        Llg.step(state);
        time+=af::timer::stop(t2);
        calcm(state, stream);
    }
    std::cout<<"Llg.get_time_heff       [af-s]: "<< Llg.get_time_heff() <<std::endl;
    std::cout<<"time                    [af-s]: "<< time <<std::endl;
    std::cout<<"Time integrate 1ns      [af-s]: "<< af::timer::stop(t) <<std::endl;
    std::cout<<"Llg.get_timer_allsteps  [af-s]: "<< Llg.get_time_allsteps() <<std::endl;
    vti_writer_micro(state.m, mesh , (filepath + "2ns").c_str());
    stream.close();
    return 0;
}

void calcm(State state, std::ostream& myfile){
    array mean_dim3 = mean(mean(mean(state.m, 0), 1), 2);
    myfile << std::setw(12) << state.t << "\t" << afvalue(mean_dim3(span, span, span, 0))<< "\t" << afvalue(mean_dim3(span, span, span, 1)) << "\t" << afvalue(mean_dim3(span, span, span, 2))<< std::endl;
}
