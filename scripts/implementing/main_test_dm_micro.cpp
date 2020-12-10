#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace magnumafcpp;

using namespace af;

int main(int argc, char** argv) {
    // const double x = 1.e-0, y = 1.e-0, z = 1.e-0;
    const int nx = 4, ny = 4, nz = 4;
    Mesh mesh(nx, ny, nz, 1, 1, 1);
    Material material = Material();
    state.Ms = 1;
    material.D = 1;
    material.D_axis[0] = 1;
    material.D_axis[1] = 0;
    material.D_axis[2] = 0;
    material.A = 1;

    // array m = constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, f64);
    // m(span, span, span, 0) = constant(1.0, mesh.nx, mesh.ny, mesh.nz, 1,
    // f64);
    std::cout << mesh.nx << std::endl;
    array m = constant(0, dim4(mesh.nx, mesh.ny, mesh.nz, 3), f64);
    m(span, span, span, 0) = 1.;
    // array m = iota(dim4(mesh.nx, mesh.ny, mesh.nz, 3), dim4(1, 1, 1, 1),
    // f64);
    print("m", m);
    // m=reorder(m, 2, 1, 0, 3);
    State state(mesh, material, m);
    std::vector<uptr_FieldTerm> llgterm;
    llgterm.push_back(uptr_FieldTerm(new DmiField(mesh, material)));
    LLG llg(state, llgterm);

    print("DmiField", llg.Fieldterms[0]->h(state));
    return 0;
}
