#include "magnum_af.hpp"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

using namespace magnumafcpp;

using namespace af;
typedef std::shared_ptr<LLGTerm> llgt_ptr;
int main(int argc, char** argv) {
    std::cout << "argc" << argc << std::endl;
    for (int i = 0; i < argc; i++)
        cout << "Parameter " << i << " was " << argv[i] << "\n";

    std::string filepath(argc > 1 ? argv[1] : "../Data/Testing");
    if (argc > 0)
        filepath.append("/");
    std::cout << "Writing into path " << filepath.c_str() << std::endl;

    setDevice(argc > 2 ? std::stoi(argv[2]) : 0);
    info();

    // Parameter initialization
    const int nx = 5, ny = 5, nz = 1; // nz=5 -> lz=(5-1)*dx
    const double dx = 1.e-10;
    // Generating Objects
    Mesh mesh(nx, ny, nz, dx, dx, dx);
    Material material = Material();
    material.p = 1.;
    // material.alpha = 1.;
    material.J_atom = 1;

    // Initial magnetic field
    array m = constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, f64);
    // m(span, span, span, 2) = -1;
    m(1, 1, 0) = 1;

    State state(mesh, material, m);
    vti_writer_atom(state.m, mesh, (filepath + "minit").c_str());

    std::vector<llgt_ptr> llgterm;

    llgterm.push_back(llgt_ptr(new AtomisticExchangeField(mesh)));
    af::print("llgterm", llgterm[0]->h(state));
    // LLG Llg(state, llgterm);

    return 0;
}
