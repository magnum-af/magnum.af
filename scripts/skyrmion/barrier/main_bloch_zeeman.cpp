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
    // if(argc>1) setDevice(std::stoi(argv[2]));
    info();

    // Parameter initialization
    const int nx = 30, ny = 30, nz = 1;
    const double dx = 1e-9;

    double n_interp = 60;
    double string_dt = 1e-13;

    // Generating Objects
    Mesh mesh(nx, ny, nz, dx, dx, dx);
    Material material = Material();
    state.Ms = 1.1e6;
    material.A = 1.6e-11;
    material.alpha = 1;
    material.D = 2 * 5.76e-3;
    material.Ku1 = 6.4e6;

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

    array zee = constant(0.0, 1, 1, 1, 3, f64);
    zee(0, 0, 0, 2) = (argc > 3 ? std::stod(argv[3]) / constants::mu0 : 0. / constants::mu0);
    af::print("zee_pre_tile", zee);
    zee = tile(zee, mesh.nx, mesh.ny, mesh.nz);

    LLGIntegrator Llg;
    // demag?//llgterm.push_back( llgt_ptr (new
    // AtomisticDipoleDipoleField(mesh)));
    Llg.llgterms.push_back(uptr_Fieldterm(new AtomisticExchangeField(mesh)));
    Llg.llgterms.push_back(uptr_Fieldterm(new AtomisticDmiField(mesh, material)));
    Llg.llgterms.push_back(uptr_Fieldterm(new AtomisticUniaxialAnisotropyField(mesh, material)));
    Llg.llgterms.push_back(uptr_Fieldterm(new ExternalField(zee)));

    Llg.relax(state);
    vti_writer_micro(state.m, mesh, filepath + "relax");
    state.t = 0;

    array last = constant(0, dims_vector(mesh), f64);
    last(span, span, span, 2) = 1;

    std::vector<State> inputimages;
    inputimages.push_back(state);
    inputimages.push_back(State(mesh, material, last));

    StringMethod string(state, inputimages, n_interp, string_dt, Llg.llgterms);
    string.run(filepath);
    return 0;
}
