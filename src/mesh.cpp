#include "mesh.hpp"

Mesh::Mesh (int inn0, int inn1, int inn2, double indx, double indy, double indz):
             n0(inn0), n1(inn1), n2(inn2),    dx(indx),    dy(indy),    dz(indz), 
             n0_exp(2*n0), n1_exp(2*n1), n2_exp((n2 == 1)? 1 : 2*n2)
{
  V = dx * dy * dz;
  dims=af::dim4(n0,n1,n2,3);
  dims_expanded=af::dim4(n0_exp,n1_exp,n2_exp,3);
}

void Mesh::print(std::ostream& stream){
    stream << "n0=" << n0 << " n1=" << n1 << " n2=" << n2 << " dx=" << dx << " dy=" << dy << " dz=" << dz << " V=" << V << " n0_exp=" << n0_exp << " n1_exp=" << n1_exp << " n2_exp=" << n2_exp << std::endl;
}
