#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace magnumafcpp;

using namespace af;


int main(int argc, char** argv) {

    std::cout << "argc = " << argc << std::endl;
    for (int i = 0; i < argc; i++)
        cout << "Parameter " << i << " was " << argv[i] << "\n";

    std::string filepath(argc > 1 ? argv[1] : "../Data/skyrmion_stoch");
    if (argc > 0)
        filepath.append("/");
    std::cout << "Writing into path " << filepath.c_str() << std::endl;

    setDevice(argc > 2 ? std::stoi(argv[2]) : 0);
    info();

    // Parameter initialization
    double length = 90e-9; //[nm]
    const double dx = 0.5e-9;
    const int nx = (int)(length / dx);
    std::cout << "nx = " << nx << std::endl;

    // Generating Objects
    Mesh mesh(nx, nx, 1, dx, dx, dx);
    Material material = Material();
    state.Ms = 580000;
    material.A = 15e-12;
    material.alpha = 1;
    material.D = 3e-3;
    material.Ku1 = 0.6e6;

    material.J_atom = 2. * material.A * dx;
    material.D_atom = material.D * pow(dx, 2);
    material.K_atom = material.Ku1 * pow(dx, 3);
    material.p = state.Ms * pow(dx, 3); // Compensate nz=1 instead of nz=4

    // Initial magnetic field
    array m = constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, f64);
    m(span, span, span, 2) = -1;
    for (int ix = 0; ix < mesh.nx; ix++) {
        for (int iy = 0; iy < mesh.ny; iy++) {
            const double rx = double(ix) - mesh.nx / 2.;
            const double ry = double(iy) - mesh.ny / 2.;
            const double r = sqrt(pow(rx, 2) + pow(ry, 2));
            if (r > nx / 4.)
                m(ix, iy, span, 2) = 1.;
        }
    }

    State state(mesh, material, m);
    vti_writer_atom(state.m, mesh, (filepath + "minit").c_str());

    std::vector<uptr_FieldTerm> llgterm;
    //    llgterm.push_back( uptr_FieldTerm (new AtomisticDipoleDipoleField(mesh)));
    llgterm.push_back(uptr_FieldTerm(new AtomisticExchangeField(mesh)));
    llgterm.push_back(uptr_FieldTerm(new AtomisticDmiField(mesh, material)));
    llgterm.push_back(uptr_FieldTerm(new AtomisticUniaxialAnisotropyField(mesh, material)));

    LLG Llg(state, llgterm);

    std::cout << "mrelax.vti not found, starting relaxation" << std::endl;
    timer t = af::timer::start();
    while (state.t < 15.e-10) {
        state.m = Llg.step(state);
    }
    double timerelax = af::timer::stop(t);
    vti_writer_atom(state.m, mesh, filepath + "relax");

    std::cout << "timerelax [af-s]: " << timerelax << " for " << Llg.counter_accepted + Llg.counter_reject
              << " steps, thereof " << Llg.counter_accepted << " Steps accepted, " << Llg.counter_reject
              << " Steps rejected" << std::endl;
    return 0;
}
