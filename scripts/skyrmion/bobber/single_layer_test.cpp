#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace magnumafcpp;

using namespace af;

int main(int argc, char** argv) {

    std::cout << "argc = " << argc << std::endl;
    for (int i = 0; i < argc; i++)
        std::cout << "Parameter " << i << " was " << argv[i] << "\n";

    std::string filepath(argc > 1 ? argv[1] : "./run/");
    if (argc > 1)
        filepath.append("/");
    std::cout << "Writing into path " << filepath.c_str() << std::endl;

    setDevice(argc > 2 ? std::stoi(argv[2]) : 0);
    info();

    // Parameter initialization
    const double x = 400e-9;
    const double y = 400e-9;
    const double z = 3e-9;

    const int nx = 128, ny = 128, nz = 1;
    const double dx = x / nx;
    const double dy = y / ny;
    const double dz = z / nz;

    // SK layer params
    const double Ms = 1371e3;  // A/m
    const double A = 15e-12;   // J/m
    const double Ku = 1.411e6; // J/m^3
    const double D = 2.5e-3;   // J/m^2
    const double Hz = 130e-3 / constants::mu0;

    // Generating Objects
    Mesh mesh(nx, ny, nz, dx, dy, dz);

    // Initial magnetic field
    array m = constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, f64);
    m(af::span, af::span, af::span, 2) = -1;
    for (int ix = 0; ix < mesh.nx; ix++) {
        for (int iy = 0; iy < mesh.ny; iy++) {
            const double rx = double(ix) - mesh.nx / 2.;
            const double ry = double(iy) - mesh.ny / 2.;
            const double r = sqrt(pow(rx, 2) + pow(ry, 2));
            if (r > nx / 4.)
                m(ix, iy, af::span, 2) = 1.;
        }
    }

    State state(mesh, Ms, m);
    state.write_vti(filepath + "minit");

    // defining interactions
    auto demag = uptr_FieldTerm(new DemagField(mesh, true, true, 0));
    auto exch = uptr_FieldTerm(new ExchangeField(A));
    auto aniso = uptr_FieldTerm(new UniaxialAnisotropyField(Ku, {0, 0, 1}));

    Material material = Material();
    material.D = D;
    material.D_axis[2] = -1;
    auto dmi = uptr_FieldTerm(new DmiField(mesh, material));

    array zee = constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, f64);
    zee(af::span, af::span, af::span, 2) = Hz;
    auto external = uptr_FieldTerm(new ExternalField(zee));

    af::print("dmi", dmi->h(state));

    LLGIntegrator Llg(1, {demag, exch, aniso, dmi, external});
    // LLGIntegrator Llg(1, {demag, exch, aniso, dmi, external});
    while (state.t < 3e-9) {
        if (state.steps % 100 == 0)
            state.write_vti(filepath + "m_step" + std::to_string(state.steps));
        Llg.step(state);
        std::cout << state.steps << "\t" << state.t << "\t" << state.meani(2) << "\t" << Llg.E(state) << std::endl;
    }
    //    Llg.relax(state);
    state.write_vti(filepath + "m_relaxed");

    // preparing string method
    //    double n_interp = 60;
    //    double string_dt=1e-13;
    //    const int string_steps = 10000;

    //    array last   = constant( 0, dims_vector(mesh), f64);
    //    last(span, span, span, 2)=1;
    //
    //    std::vector<State> inputimages;
    //    inputimages.push_back(state);
    //    inputimages.push_back(State(mesh, material, last));
    //
    //    StringMethod string(state, inputimages, n_interp, string_dt , Llg.llgterms);
    //    string.run(filepath);
    return 0;
}
