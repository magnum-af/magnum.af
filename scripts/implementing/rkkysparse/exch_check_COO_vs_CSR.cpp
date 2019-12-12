#include "magnum_af.hpp"
#include "arrayfire.h"
#include <iostream>

using namespace magnumafcpp;
using namespace af;

int main(int argc, char **argv)
{

    std::cout << "argc = " << argc << std::endl;
    for (int i = 0; i < argc; i++)
        std::cout << "Parameter " << i << " was " << argv[i] << "\n";
    std::string filepath(argc > 1 ? argv[1] : "./run/");
    if (argc > 1)
        filepath.append("/");
    std::cout << "Writing into path " << filepath.c_str() << std::endl;
    setDevice(argc > 2 ? std::stoi(argv[2]) : 0);
    info();

    const int nx = 30, ny = 40, nz = 5;
    const double dx = 1e-9, dy = 2e-9, dz = 3e-9;
    array A = af::randu(nx, ny, nz, 3, f64);
    A = A * 1e-12;
    //Note: uncomment for const A check//const double A = 2.3e-12;
    Mesh mesh(nx, ny, nz, dx, dy, dz);
    auto coo = SparseExchangeField(A, mesh, true, true);
    af::array coo_to_dense = af::sparseConvertTo(coo.matr, AF_STORAGE_DENSE);
    auto csr = SparseExchangeField(A, mesh, true, false);
    af::array csr_to_dense = af::sparseConvertTo(csr.matr, AF_STORAGE_DENSE);

    af::print("logic max", af::max(af::flat(csr_to_dense == coo_to_dense)));
    af::print("logic min", af::min(af::flat(csr_to_dense == coo_to_dense)));

    af::print("max  csr", af::max(af::flat(csr_to_dense)));
    af::print("mean csr", af::mean(af::flat(csr_to_dense)));
    af::print("min  csr", af::min(af::flat(csr_to_dense)));

    af::array csr_max_val;
    af::array csr_max_idx;
    af::max(csr_max_val, csr_max_idx, af::flat(csr_to_dense));
    af::print("csr_max_val", csr_max_val);
    af::print("csr_max_idx", csr_max_idx);

    af::print("max  coo", af::max(af::flat(coo_to_dense)));
    af::print("mean coo", af::mean(af::flat(coo_to_dense)));
    af::print("min  coo", af::min(af::flat(coo_to_dense)));

    af::array coo_max_val;
    af::array coo_max_idx;
    af::max(coo_max_val, coo_max_idx, af::flat(coo_to_dense));
    af::print("coo_max_val", coo_max_val);
    af::print("coo_max_idx", coo_max_idx);

    af::print("matr max  diff", af::max(af::flat(csr_to_dense - coo_to_dense)));
    af::print("matr mean diff", af::mean(af::flat(csr_to_dense - coo_to_dense)));
    af::print("matr min  diff", af::min(af::flat(csr_to_dense - coo_to_dense)));
    if (afvalue(csr_max_val) == afvalue(coo_max_val))
        std::cout << " maxvals equal" << std::endl;
    else
        std::cout << " maxvals not equal" << std::endl;
    if (afvalue(af::mean(af::flat(coo_to_dense))) == afvalue(af::mean(af::flat(csr_to_dense))))
        std::cout << " mean equal" << std::endl;
    else
        std::cout << " mean not equal" << std::endl;

    af::array m = af::randu(nx, ny, nz, 3, f64);
    State state(mesh, 1e6, m);
    af::array hcsr = csr.h(state);
    af::array hcoo = coo.h(state);
    af::print("1e20 * max  diff h", 1e20 * af::max(af::flat(hcsr - hcoo)));
    af::print("1e20 * mean diff h", 1e20 * af::mean(af::flat(hcsr - hcoo)));
    af::print("1e20 * min  diff h", 1e20 * af::min(af::flat(hcsr - hcoo)));
    std::cout.precision(12);
    std::cout << "Ecoo=" << coo.E(state) << std::endl;
    std::cout << "Ecsr=" << csr.E(state) << std::endl;
    std::cout << "Ecoo=" << coo.E(state, hcoo) << std::endl;
    std::cout << "Ecsr=" << csr.E(state, hcsr) << std::endl;

    return 0;
}
