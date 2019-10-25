#include "magnum_af.hpp"
#include "arrayfire.h"
#include <iostream>

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


    //const int nx=4, ny=4 , nz=5;
    //af::array values = af::iota(af::dim4(nx * ny * nz * 3, 1), af::dim4(1), f64);
    //af::print("values", values);
    //af::array val2 = af::moddims(values, nx, ny, nz, 3);
    //af::print("val2", val2);
    //af::print("", val2(1,0,1,1));
    //af::print("", val2(1,1,1,1));
    //af::print("", val2(1,2,1,1));

    //const int nx=4, ny=4 , nz=5;
    //af::array A = af::constant(0, nx, nz, ny, 3, f64);
    //A(0,1,0,0)=1.;
    //A(0,1,2,0)=1.;
    //A(3,1,0,0)=1.;
    //A(3,0,1,0)=1.;
    //af::print("", A);

    //const int nx=4, ny=4;
    //af::array A = af::constant(0, nx, ny, f64);
    //A(0,1)=1.;
    //A(1,0)=2.;
    //A(2,2)=3.;
    //A(3,0)=4.;
    //A(3,3)=5.;

    //af::print("", A);
    //af::array Asparse = af::sparse(A);
    //af::print("", Asparse);

    //af::array vals = af::iota(af::dim4(5), af::dim4(1), f64);
    //vals = vals + 1;

    //// Moving arraz to COO to CSR to dense
    //std::vector<double> cvals   = {1,2,3,4,5,6,7};
    //std::vector<int> cnindex = {0,1,2,3,4,5,6};
    //std::vector<int> cmindex = {0,1,2,3,4,5,6};
    //af::array vals = af::array(cvals.size(), cvals.data());
    //af::array mindex = af::array(cmindex.size(), cmindex.data());
    //af::array nindex = af::array(cnindex.size(), cnindex.data());
    //af::print("vals", vals);
    //af::print("vals", mindex);
    //af::print("vals", nindex);
    //af::array A_COO = af::sparse((dim_t) cvals.size(), (dim_t) cvals.size(), vals, mindex, nindex, AF_STORAGE_COO);
    //af::print("", A_COO);
    //af::array A_CSR = af::sparseConvertTo(A_COO, AF_STORAGE_CSR);
    //af::print("", A_CSR);
    //af::array A_dense = af::sparseConvertTo(A_CSR, AF_STORAGE_DENSE);
    //af::print("", A_dense);

    //const int nx =20, ny=30 , nz=40;
    //max diff nonzero// const int nx = 30, ny=4 , nz=5;
    //2nd field crashes with invalid buffer size//const int nx = 128, ny=128 , nz=4;
    const int nx = 30, ny=4 , nz=5;
    //const int nx = 30, ny=40 , nz=5;
    //good numerics//const int nx = 30, ny=11 , nz=12;
    const double dx = 1e-9, dy=2e-9, dz=3e-9;
    //array A = af::constant(1e-5, nx, ny, nz, 3, f64);
    array A = af::randu(nx, ny, nz, 3, f64);
    //array A = af::iota(af::dim4(nx, ny, nz, 3), af::dim4(1,1,1,1), f64);
    //A = A + 1;
    A = A * 1e-12;//(nx*nz*ny);
    array RKKY = af::constant(0, nx, ny, nz, 3, f64);
    Mesh mesh(nx, ny, nz, dx, dy, dz);
    auto coo = RKKYExchangeField(int(2), RKKY_values(RKKY), Exchange_values(A), mesh, af::array());
    af::array coo_to_dense = af::sparseConvertTo(coo.matr, AF_STORAGE_DENSE);
    auto csr = RKKYExchangeField(RKKY_values(RKKY), Exchange_values(A), mesh);
    af::array csr_to_dense = af::sparseConvertTo(csr.matr, AF_STORAGE_DENSE);
    //af::print("", exch.matr);
    //af::print("", csr_to_dense);
    //af::print("", coo_to_dense);
    //af::print("", csr_to_dense == coo_to_dense);

    af::print("logic max", af::max(af::flat(csr_to_dense == coo_to_dense)));
    af::print("logic min", af::min(af::flat(csr_to_dense == coo_to_dense)));

    af::print("max  csr", af::max (af::flat(csr_to_dense)));
    af::print("mean csr", af::mean(af::flat(csr_to_dense)));
    af::print("min  csr", af::min (af::flat(csr_to_dense)));

    af::array csr_max_val;
    af::array csr_max_idx;
    af::max(csr_max_val, csr_max_idx, af::flat(csr_to_dense));
    af::print("csr_max_val", csr_max_val);
    af::print("csr_max_idx", csr_max_idx);

    af::print("max  coo", af::max (af::flat(coo_to_dense)));
    af::print("mean coo", af::mean(af::flat(coo_to_dense)));
    af::print("min  coo", af::min (af::flat(coo_to_dense)));

    af::array coo_max_val;
    af::array coo_max_idx;
    af::max(coo_max_val, coo_max_idx, af::flat(coo_to_dense));
    af::print("coo_max_val", coo_max_val);
    af::print("coo_max_idx", coo_max_idx);

    af::print("max  diff", af::max (af::flat(csr_to_dense - coo_to_dense)));
    af::print("mean diff", af::mean(af::flat(csr_to_dense - coo_to_dense)));
    af::print("min  diff", af::min (af::flat(csr_to_dense - coo_to_dense)));
    //if ( csr == coo ) std::cout << "csr == coo" << std::endl;
    //if ( csr != coo ) std::cout << "csr != coo" << std::endl;
    if (afvalue(csr_max_val) == afvalue(coo_max_val)) std::cout << " maxvals equal" << std::endl; else std::cout << " maxvals not equal" << std::endl;
    if (afvalue(af::mean(af::flat(coo_to_dense))) == afvalue(af::mean(af::flat(csr_to_dense)))) std::cout << " mean equal" << std::endl; else std::cout << " mean not equal" << std::endl;

    af::array m = af::randu(nx, ny, nz, 3, f64);
    State state (mesh, 1e6, m);
    af::array hcsr = csr.h(state);
    af::array hcoo = coo.h(state);
    af::print("max  diff", 1e20 * af::max (af::flat(hcsr - hcoo)));
    af::print("mean diff", 1e20 * af::mean(af::flat(hcsr - hcoo)));
    af::print("min  diff", 1e20 * af::min (af::flat(hcsr - hcoo)));

    return 0;
}
