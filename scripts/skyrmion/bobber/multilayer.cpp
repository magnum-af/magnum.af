#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace magnumafcpp;

using namespace af;

int main(int argc, char** argv)
{

    std::cout<<"argc = "<<argc<<std::endl;
     for (int i=0; i<argc; i++)
          std::cout << "Parameter " << i << " was " << argv[i] << "\n";

    std::string filepath(argc>1? argv[1]: "./run/");
    if(argc>1)filepath.append("/");
    std::cout<<"Writing into path "<<filepath.c_str()<<std::endl;

    setDevice(argc>2? std::stoi(argv[2]):0);
    info();

    // Parameter initialization
    const double x=400e-9;
    const double y=400e-9;
    const double z=34e-9;
    //TODO fix RKKY adjacent layers problem and move from 5 4 layers

    //const int nx = 1, ny=1 , nz=34;//TODO
    const int nx = 100, ny=100 , nz=34;//TODO this discretiyation stabiliyes skym with scaled RKKYval=0.8e-3 * 1e-9
    //const int nx = 128, ny=128 , nz=4;
    const double dx= x/nx;
    const double dy= y/ny;
    const double dz= z/nz;

    const double Hz = 130e-3/constants::mu0;
    const double RKKY_val = 0.8e-3 * 1e-9;
    //TODO maybe 0.5 factor (see mumax3):
    //const double RKKY_val = 0.8e-3 * 1e-9* 0.5;//
    //NOTE//const double RKKY_val = 0.8e-3;//NOTE: causes NaNs during integration
    // SK layer params
    const double SK_Ms =1371e3;// A/m
    const double SK_A = 15e-12;// J/m
    const double SK_Ku = 1.411e6;// J/m^3
    const double SK_D =2.5e-3;// J/m^2

    // Ferrimagnetic interface layer params
    const double IL_Ms = 488.2e3;// A/m
    const double IL_A = 4e-12;// J/m
    const double IL_Ku = 486.6e3;// J/m^3

    array geom = af::constant(0.0, nx, ny, nz, 3, f64);
    geom(af::span, af::span,  0, af::span) = 1;
    geom(af::span, af::span,  3, af::span) = 1;
    geom(af::span, af::span,  6, af::span) = 1;
    geom(af::span, af::span,  9, af::span) = 1;
    geom(af::span, af::span, 12, af::span) = 1;

    geom(af::span, af::span, 15, af::span) = 2;
    geom(af::span, af::span, 16, af::span) = 2;
    geom(af::span, af::span, 17, af::span) = 2;
    geom(af::span, af::span, 18, af::span) = 2;

    geom(af::span, af::span, 19, af::span) = 1;
    geom(af::span, af::span, 22, af::span) = 1;
    geom(af::span, af::span, 25, af::span) = 1;
    geom(af::span, af::span, 28, af::span) = 1;
    geom(af::span, af::span, 31, af::span) = 1;

    //af::print("geom", geom);
    //af::print("geom == 1", geom == 1);
    //af::print("geom == 2", geom == 2);

    af::array todouble = af::constant(1., nx, ny, nz, 3, f64);
    af::array Ms = ( SK_Ms * (geom == 1) + IL_Ms * (geom == 2)) * todouble;
    af::array A  = ( SK_A  * (geom == 1) + IL_A  * (geom == 2)) * todouble;
    af::array Ku = ( SK_Ku * (geom == 1) + IL_Ku * (geom == 2)) * todouble;
    af::array D  = ( SK_D  * (geom == 1)) * todouble; // + IL_D * (geom == 2);
    std::cout << "type" << geom.type() << " " << Ms.type() << std::endl;

    //af::print("Ms", Ms);
    //af::print("1e12 * A", 1e12 * A);
    //af::print("Ku", Ku);
    //af::print("D", D);

    array RKKY = af::constant(0.0, nx, ny, nz, 3, f64);
    RKKY(af::span, af::span, 18, af::span) = RKKY_val;
    RKKY(af::span, af::span, 19, af::span) = RKKY_val;

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
    m(af::span, af::span, 1, af::span) = 0;
    m(af::span, af::span, 2, af::span) = 0;

    m(af::span, af::span, 4, af::span) = 0;
    m(af::span, af::span, 5, af::span) = 0;

    m(af::span, af::span, 7, af::span) = 0;
    m(af::span, af::span, 8, af::span) = 0;

    m(af::span, af::span, 10, af::span) = 0;
    m(af::span, af::span, 11, af::span) = 0;

    m(af::span, af::span, 13, af::span) = 0;
    m(af::span, af::span, 14, af::span) = 0;

    m(af::span, af::span, 20, af::span) = 0;
    m(af::span, af::span, 21, af::span) = 0;

    m(af::span, af::span, 23, af::span) = 0;
    m(af::span, af::span, 24, af::span) = 0;

    m(af::span, af::span, 26, af::span) = 0;
    m(af::span, af::span, 27, af::span) = 0;

    m(af::span, af::span, 29, af::span) = 0;
    m(af::span, af::span, 30, af::span) = 0;

    m(af::span, af::span, 32, af::span) = 0;
    m(af::span, af::span, 33, af::span) = 0;

    State state(mesh, Ms, m);
    state.write_vti(filepath + "minit");

    // defining interactions
    auto demag = LlgTerm (new DemagField(mesh, true, true, 0));
    auto exch = LlgTerm (new RKKYExchangeField(RKKY_values(RKKY), Exchange_values(A), mesh));
    auto aniso = LlgTerm (new UniaxialAnisotropyField(Ku, (std::array<double ,3>) {0, 0, 1}));

    //NOTE try//auto dmi = LlgTerm (new DmiField(SK_D, {0, 0, -1}));
    auto dmi = LlgTerm (new DmiField(D, {0, 0, -1}));

    array zee = constant(0.0, mesh.n0, mesh.n1, mesh.n2, 3, f64);
    zee(af::span, af::span, af::span, 2) = Hz;
    auto external = LlgTerm (new ExternalField(zee));

    //af::print("dmi", dmi->h(state));
    //af::print("exch", exch->h(state));

    LLGIntegrator Llg(1, {demag, exch, aniso, dmi, external});
    //LLGIntegrator Llg(1, {demag, exch, aniso, dmi, external});
    while (state.t < 3e-9){
        if (state.steps % 100 == 0) state.write_vti(filepath + "m_step" + std::to_string(state.steps));
        Llg.step(state);
        std::cout << state.steps << "\t" << state.t << "\t" <<  state.meani(2) << "\t" << Llg.E(state) << std::endl;
    }
    //Llg.relax(state);
    state.write_vti(filepath + "m_relaxed");

//  // preparing string method
//    double n_interp = 60;
//    double string_dt=1e-13;
//    const int string_steps = 10000;
//
//    array last   = constant( 0, mesh.dims, f64);
//    last(span, span, span, 2)=1;
//
//    std::vector<State> inputimages;
//    inputimages.push_back(state);
//    inputimages.push_back(State(mesh, material, last));
//
//    String string(state, inputimages, n_interp, string_dt , Llg.llgterms);
//    string.run(filepath);
    return 0;
}
