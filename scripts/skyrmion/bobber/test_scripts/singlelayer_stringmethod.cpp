#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace magnumafcpp;

using namespace af;

int main(int argc, char** argv)
{

    std::cout<<"argc = "<<argc<<std::endl;
    for (int i=0; i<argc; i++) std::cout << "Parameter " << i << " was " << argv[i] << "\n";
    std::string filepath(argc>1? argv[1]: "./run/");
    if(argc>1)filepath.append("/");
    std::cout<<"Writing into path "<<filepath.c_str()<<std::endl;

    setDevice(argc>2? std::stoi(argv[2]):0);
    info();
    StageTimer timer;

    // Parameter initialization
    const int nx = 100, ny=100 , nz=1;
    const double x=400e-9;
    const double y=400e-9;
    const double z= nz * 1e-9;

    const double dx= x/nx;
    const double dy= y/ny;
    const double dz= z/nz;

    const double Hz = 130e-3/constants::mu0;
    //const double RKKY_val = 0.8e-3 * 1e-9* 0.5;

    // SK layer params
    const double Ms =1371e3;// A/m
    const double A = 15e-12;// J/m
    const double Ku = 1.411e6;// J/m^3
    const double D =2.5e-3;// J/m^2

    //Generating Objects
    Mesh mesh(nx, ny, nz, dx, dy, dz);

    // Initial magnetic field
    array m = constant(0.0, mesh.n0, mesh.n1, mesh.n2, 3, f64);
    m(af::span, af::span, af::span, 2) = -1;
    for(int ix=0;ix<mesh.n0;ix++){
        for(int iy=0;iy<mesh.n1;iy++){
            const double rx=double(ix)-mesh.n0/2.;
            const double ry=double(iy)-mesh.n1/2.;
            const double r = sqrt(pow(rx, 2)+pow(ry, 2));
            if(r>nx/4.) m(ix, iy, af::span, 2)=1.;
        }
    }

    // defining interactions
    auto demag = LlgTerm (new DemagField(mesh, true, true, 0));
    auto exch = LlgTerm (new ExchangeField(A));
    auto aniso = LlgTerm (new UniaxialAnisotropyField(Ku, (std::array<double ,3>) {0, 0, 1}));
    auto dmi = LlgTerm (new DmiField(D, {0, 0, -1}));
    af::array zee = af::constant(0.0, mesh.n0, mesh.n1, mesh.n2, 3, f64);
    zee(af::span, af::span, af::span, 2) = Hz;
    auto external = LlgTerm (new ExternalField(zee));

    //af::print("dmi", dmi->h(state_1));
    //af::print("exch", exch->h(state_1));

    LLGIntegrator llg(1, {demag, exch, aniso, dmi, external});
    timer.print_stage("init ");

    State state_1(mesh, Ms, m);

    if (! exists(filepath + "m_relaxed.vti" ) ){
        std::cout << "Relaxing minit" << std::endl;
        state_1.write_vti(filepath + "minit");
        //LLGIntegrator Llg(1, {demag, exch, aniso, dmi, external});
        while (state_1.t < 3e-9){
            if (state_1.steps % 100 == 0) state_1.write_vti(filepath + "m_step" + std::to_string(state_1.steps));
            llg.step(state_1);
            std::cout << std::scientific << state_1.steps << "\t" << state_1.t << "\t" <<  state_1.meani(2) << "\t" << llg.E(state_1) << std::endl;
        }
        //Llg.relax(state_1);
        timer.print_stage("relax");
        state_1.write_vti(filepath + "m_relaxed");
    }
    else {
        std::cout << "Found m_relaxed, reading in state_1." << std::endl;
        state_1._vti_reader(filepath + "m_relaxed.vti");
        state_1.write_vti(filepath + "m_relaxed_from_read_in");

    }

    // preparing string method
    double n_interp = 60;
    double string_dt=1e-13;
    //const int string_steps = 1000;

    af::array m2 = af::constant(0.0, mesh.n0, mesh.n1, mesh.n2, 3, f64);
    m2(af::span, af::span, af::span, 2) = 1;

    State state_2(mesh, Ms, m2);

    vti_writer_micro(state_2.m, state_2.mesh, filepath + "m_state_2_init");

    while (state_2.t < 3e-9){
        if (state_2.steps % 100 == 0) state_2.write_vti(filepath + "m_state2_relaxing" + std::to_string(state_2.steps));
        llg.step(state_2);
        std::cout << std::scientific << state_2.steps << "\t" << state_2.t << "\t" <<  state_2.meani(2) << "\t" << llg.E(state_2) << std::endl;
    }
    vti_writer_micro(state_2.m, state_2.mesh, filepath + "m_state_2_relaxed");
    timer.print_stage("state2 relaxed");

    std::vector<State> inputimages;
    inputimages.push_back(state_1);
    inputimages.push_back(state_2);

    String string(state_1, inputimages, n_interp, string_dt , llg);
    string.run(filepath);
    timer.print_stage("string relaxed");
    return 0;
}
