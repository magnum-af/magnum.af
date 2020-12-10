#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace magnumafcpp;

using namespace af;

bool compare(double a, double b) {
    if (fabs(a - b) / fabs(a + b) < 1e-30)
        return false;
    else
        return true;
}

int main(int argc, char** argv) {
    info();
    std::string filepath(argc > 0 ? argv[1] : "../Data/Testing/");
    std::cout << filepath << std::endl;
    int nx = 2, ny = 1, nz = 1; // nz=5 -> lz=(5-1)*dx
    // const double dx=1;
    const double dx = 2.715e-10;

    // Generating Objects
    Mesh mesh(nx, ny, nz, dx, dx, dx);
    Material material = Material();
    material.alpha = 1;
    state.material.afsync = false;
    // material.p=1;
    material.p = 9.274009994e-24;

    // Initial magnetic field
    //
    //-------------------------------------------------------
    array m = constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, f64);
    m(0, 0, 0, 0) = 0;
    m(0, 0, 0, 1) = 0;
    m(0, 0, 0, 2) = 1;

    m(1, 0, 0, 0) = 0;
    m(1, 0, 0, 1) = 0;
    m(1, 0, 0, 2) = 1;
    State state(mesh, material, m);
    vti_writer_atom(state.m, mesh, (filepath + "/minit").c_str());

    std::vector<uptr_FieldTerm> llgterm;
    llgterm.push_back(uptr_FieldTerm(new AtomisticDipoleDipoleField(mesh)));
    LLG llg(state, llgterm);
    double analytical = -pow(material.p, 2) * constants::mu0 / (4. * M_PI) / pow(dx, 3);
    // std::cout << "ENERGY    = " << llg.E(state) <<std::endl;
    std::cout << "Analytical= " << analytical << std::endl;
    if (compare(llg.E(state), analytical))
        std::cout << "!!! TEST FAILED !!!" << std::endl;
    // std::cout << "Analytical= " << - pow(material.p,
    // 2)*constants::mu0/(4.*M_PI)/pow(dx, 3) <<std::endl;
    std::cout << "H_dip_1   = " << 0 << ", " << 0 << ", " << material.p / (4 * M_PI * pow(dx, 3)) << std::endl;
    std::cout << "H_dip_2   = " << 0 << ", " << 0 << ", " << material.p / (4 * M_PI * pow(dx, 3)) << std::endl;

    //-------------------------------------------------------
    m = constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, f64);
    m(0, 0, 0, 0) = 0;
    m(0, 0, 0, 1) = 0;
    m(0, 0, 0, 2) = 1;

    m(1, 0, 0, 0) = 1;
    m(1, 0, 0, 1) = 0;
    m(1, 0, 0, 2) = 0;

    state.m = m;
    vti_writer_atom(state.m, mesh, (filepath + "/minit2").c_str());
    std::cout << "ENERGY    = " << llg.E(state) << std::endl;
    std::cout << "Analytical= " << 0 << std::endl;
    std::cout << "H_dip_1   = " << -2 * material.p / (4 * M_PI * pow(dx, 3)) << ", " << 0 << ", " << 0 << std::endl;
    std::cout << "H_dip_2   = " << 0 << ", " << 0 << ", " << material.p / (4 * M_PI * pow(dx, 3)) << std::endl;
    //-------------------------------------------------------
    m = constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, f64);
    m(0, 0, 0, 0) = 0;
    m(0, 0, 0, 1) = 0;
    m(0, 0, 0, 2) = 1;

    m(1, 0, 0, 0) = 0;
    m(1, 0, 0, 1) = 0;
    m(1, 0, 0, 2) = -1;
    state.m = m;
    vti_writer_atom(state.m, mesh, (filepath + "/minit2").c_str());
    std::cout << "ENERGY    = " << llg.E(state) << std::endl;
    std::cout << "Analytical= " << pow(material.p, 2) * constants::mu0 / (4. * M_PI) / pow(dx, 3) << std::endl;
    std::cout << "H_dip_1   = " << 0 << ", " << 0 << ", " << -material.p / (4 * M_PI * pow(dx, 3)) << std::endl;
    std::cout << "H_dip_2   = " << 0 << ", " << 0 << ", " << material.p / (4 * M_PI * pow(dx, 3)) << std::endl;
    //-------------------------------------------------------

    nx = 1, ny = 2, nz = 1;
    mesh = Mesh(nx, ny, nz, dx, dx, dx);
    m = constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, f64);
    m(0, 0, 0, 0) = 0;
    m(0, 0, 0, 1) = 0;
    m(0, 0, 0, 2) = 1;

    m(0, 1, 0, 0) = 0;
    m(0, 1, 0, 1) = 0;
    m(0, 1, 0, 2) = 1;
    state = State(mesh, material, m);
    vti_writer_atom(state.m, mesh, (filepath + "/minit").c_str());

    llgterm.pop_back();
    llgterm.push_back(uptr_FieldTerm(new AtomisticDipoleDipoleField(mesh)));
    // TODO this leads to compiler error
    // llg=LLG(state, llgterm);
    LLG llg2(state, llgterm);
    std::cout << "ENERGY    = " << llg2.E(state) << std::endl;
    std::cout << "Analytical= " << -pow(material.p, 2) * constants::mu0 / (4. * M_PI) / pow(dx, 3)
              << std::endl; // TODO calc on paper, but should be like case 1
    //-------------------------------------------------------
    m = constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, f64);
    m(0, 0, 0, 0) = 1;
    m(0, 0, 0, 1) = 0;
    m(0, 0, 0, 2) = 0;

    m(0, 1, 0, 0) = 0;
    m(0, 1, 0, 1) = 0;
    m(0, 1, 0, 2) = 1;
    state.m = m;
    std::cout << "ENERGY    = " << llg2.E(state) << std::endl;
    std::cout << "Analytical= " << 0 << std::endl;
    return 0;
}
