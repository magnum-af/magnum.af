#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace magnumafcpp;

using namespace af;

bool compare(double a, double b) {
    // std::cout << "COM:"<< a <<", " << b <<",
    // "<<fabs(a-b)/fabs(a+b)<<std::endl;
    if (a == 0 && b == 0)
        return false;
    if (fabs(a + b) == 0) {
        std::cout << "DIVISION by 0 in compare" << std::endl;
        return true;
    }
    if (fabs(a - b) / fabs(a + b) < 1e-15)
        return false;
    else
        return true;
}

int main(int argc, char** argv) {
    info();
    std::cout.precision(32);
    int nx = 5, ny = 5, nz = 5;
    const double dx = 2.715e-10;

    // Generating Objects
    Mesh mesh(nx, ny, nz, dx, dx, dx);
    Material material = Material();
    state.Ms = 1.1e6;
    material.A = 1.6e-11;
    material.J_atom = 2. * material.A * dx;
    material.p = state.Ms * pow(dx, 3);

    //-------------------------------------------------------
    array m = randu(mesh.nx, mesh.ny, mesh.nz, 3, f64);
    State state(mesh, material, m);

    std::vector<uptr_FieldTerm> llgterm;
    llgterm.push_back(uptr_FieldTerm(new AtomisticExchangeField(mesh)));
    LLG llg(state, llgterm);
    std::vector<uptr_FieldTerm> llgterm2;
    llgterm2.push_back(uptr_FieldTerm(new ExchangeField(mesh, material)));
    LLG llg2(state, llgterm2);

    for (int x = 0; x < nx; x++) {
        for (int y = 0; y < ny; y++) {
            for (int z = 0; z < nz; z++) {
                for (int m = 0; m < 3; m++) {
                    if (compare(util::afvalue_as_f64(llg.Fieldterms[0]->H_in_Apm(state)(x, y, z, m)),
                                util::afvalue_as_f64(llg2.Fieldterms[0]->H_in_Apm(state)(x, y, z, m)))) {
                        std::cout << "!!! TEST  FAILED at " << x << " " << y << " " << z << " " << m << std::endl;
                        std::cout << util::afvalue_as_f64(llg.Fieldterms[0]->H_in_Apm(state)(x, y, z, m)) << " , "
                                  << util::afvalue_as_f64(llg2.Fieldterms[0]->H_in_Apm(state)(x, y, z, m)) << std::endl;
                    }
                }
            }
        }
    }
    return 0;
}
