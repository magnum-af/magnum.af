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
    const double z=5e-9;
    //TODO fix RKKY adjacent layers problem and move from 5 4 layers

    const int nx = 80, ny=80 , nz=5;//TODO this discretiyation stabiliyes skym with scaled RKKYval=0.8e-3 * 1e-9
    //const int nx = 128, ny=128 , nz=5;
    const double dx= x/nx;
    const double dy= y/ny;
    const double dz= z/nz;

    const double Hz = 130e-3/constants::mu0;
    //const double RKKY_val = 0.8e-3;
    const double RKKY_val = 0.8e-3 * 1e-9;//TODO
    // SK layer params
    const double SK_Ms =1371e3;// A/m
    const double SK_A = 15e-12;// J/m
    const double SK_Ku = 1.411e6;// J/m^3
    const double SK_D =2.5e-3;// J/m^2

    // Ferrimagnetic interface layer params
    const double IL_Ms = 488.2e3;// A/m
    const double IL_A = 4e-12;// J/m
    const double IL_Ku = 486.6e3;// J/m^3

    array Ms = af::constant(0.0, nx, ny, nz, 3, f64);
    Ms(af::span, af::span, 0, af::span) = SK_Ms;
    Ms(af::span, af::span, 1, af::span) = IL_Ms;
    Ms(af::span, af::span, 2, af::span) = IL_Ms;
    Ms(af::span, af::span, 3, af::span) = IL_Ms;
    Ms(af::span, af::span, 4, af::span) = SK_Ms;

    array A = af::constant(0.0, nx, ny, nz, 3, f64);
    A(af::span, af::span, 0, af::span) = SK_A;
    A(af::span, af::span, 1, af::span) = IL_A;
    A(af::span, af::span, 2, af::span) = IL_A;
    A(af::span, af::span, 3, af::span) = IL_A;
    A(af::span, af::span, 4, af::span) = SK_A;

    array Ku = af::constant(0.0, nx, ny, nz, 3, f64);
    Ku(af::span, af::span, 0, af::span) = SK_Ku;
    Ku(af::span, af::span, 1, af::span) = IL_Ku;
    Ku(af::span, af::span, 2, af::span) = IL_Ku;
    Ku(af::span, af::span, 3, af::span) = IL_Ku;
    Ku(af::span, af::span, 4, af::span) = SK_Ku;

    array D = af::constant(0.0, nx, ny, nz, 3, f64);
    D(af::span, af::span, 0, af::span) = SK_D;
    D(af::span, af::span, 4, af::span) = SK_D;


    array RKKY = af::constant(0.0, nx, ny, nz, 3, f64);
    RKKY(af::span, af::span, 0, af::span) = RKKY_val;
    RKKY(af::span, af::span, 1, af::span) = RKKY_val;
    RKKY(af::span, af::span, 2, af::span) = 0.;
    RKKY(af::span, af::span, 3, af::span) = RKKY_val;
    RKKY(af::span, af::span, 4, af::span) = RKKY_val;
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

    State state(mesh, Ms, m);
    state.write_vti(filepath + "minit");

    // defining interactions
    auto demag = LlgTerm (new DemagField(mesh, true, true, 0));
    //auto exch = LlgTerm (new ExchangeField(A));
    //TODO causes Nans//auto exch = LlgTerm (new RKKYExchangeField(RKKY_values(RKKY), Exchange_values(A), mesh));
    auto exch = LlgTerm (new RKKYExchangeField(RKKY_values(RKKY), Exchange_values(A), mesh));
    auto aniso = LlgTerm (new UniaxialAnisotropyField(Ku, (std::array<double ,3>) {0, 0, 1}));

    auto dmi = LlgTerm (new DmiField(D, {0, 0, -1}));

    array zee = constant(0.0, mesh.n0, mesh.n1, mesh.n2, 3, f64);
    zee(af::span, af::span, af::span, 2) = Hz;
    auto external = LlgTerm (new ExternalField(zee));

    //af::print("dmi", dmi->h(state));
    af::print("exch", exch->h(state));

    LLGIntegrator Llg(1, {demag, exch, aniso, dmi, external});
    //LLGIntegrator Llg(1, {demag, exch, aniso, dmi, external});
    while (state.t < 3e-9){
        if (state.steps % 100 == 0) state.write_vti(filepath + "m_step" + std::to_string(state.steps));
        Llg.step(state);
        std::cout << state.steps << "\t" << state.t << "\t" <<  state.meani(2) << "\t" << Llg.E(state) << std::endl;
    }
//    Llg.relax(state);
    state.write_vti(filepath + "m_relaxed");

    // preparing string method
//    double n_interp = 60;
//    double string_dt=1e-13;
//    const int string_steps = 10000;

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
