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
    llgterm.push_back(uptr_FieldTerm(new AtomisticDipoleDipoleField(mesh)));
    llgterm.push_back(uptr_FieldTerm(new AtomisticExchangeField(mesh)));
    llgterm.push_back(uptr_FieldTerm(new AtomisticDmiField(mesh, material)));
    llgterm.push_back(uptr_FieldTerm(new AtomisticUniaxialAnisotropyField(mesh, material)));

    array zeeswitch = constant(0.0, 1, 1, 1, 3, f64);
    zeeswitch(0, 0, 0, 2) = -0.07 * pow(state.Ms, 2) * constants::mu0;
    zeeswitch = tile(zeeswitch, mesh.nx, mesh.ny, mesh.nz);
    llgterm.push_back(uptr_FieldTerm(new ExternalField(zeeswitch)));

    LLG llg(state, llgterm);
    llg.write_fieldterms_micro(state, filepath + "init_field_micro_");
    llg.write_fieldterms_atom(state, filepath + "init_field_atom_");

    timer t = af::timer::start();
    double E_prev = 1e20;
    while (fabs((E_prev - llg.E(state)) / E_prev) > 1e-10) {
        E_prev = llg.E(state);
        for (int i = 0; i < 100; i++) {
            state.m = llg.step(state);
        }
        if (state.steps % 1000 == 0)
            std::cout << "step " << state.steps << "reldiff= " << fabs((E_prev - llg.E(state)) / E_prev) << std::endl;
    }
    double timerelax = af::timer::stop(t);
    vti_writer_atom(state.m, mesh, filepath + "relax");
    llg.write_fieldterms_micro(state, filepath);
    llg.write_fieldterms_micro(state, filepath + "field_micro_");
    llg.write_fieldterms_atom(state, filepath + "field_atom_");
    // vti_writer_atom(llg.Fieldterms[0]->h(state), mesh , filepath + "aDemag");
    // vti_writer_atom(llg.Fieldterms[1]->h(state), mesh , filepath + "aExch");
    // vti_writer_atom(llg.Fieldterms[2]->h(state), mesh , filepath + "aDMI");
    // vti_writer_atom(llg.Fieldterms[3]->h(state), mesh , filepath + "aAni");
    // vti_writer_atom(llg.Fieldterms[4]->h(state), mesh , filepath + "aZee");

    std::cout << "timerelax [af-s]: " << timerelax << " for " << llg.counter_accepted + llg.counter_reject
              << " steps, thereof " << llg.counter_accepted << " Steps accepted, " << llg.counter_reject
              << " Steps rejected" << std::endl;
    return 0;
}
