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

af::array Mesh::skyrmconf(const bool point_up){
// Returns a initial configuration to be relaxed into a skyrmion
// if point_up is true, skyrmion centers points in +z, if false in -z
     af::array m = af::constant(0.0,this->n0,this->n1,this->n2,3,f64);
     if (point_up){
         m(af::span,af::span,af::span,2) = 1.;
     }
     else {
         m(af::span,af::span,af::span,2) = -1.;
     }
     for(int ix=0;ix<this->n0;ix++){
         for(int iy=0;iy<this->n1;iy++){
             const double rx=double(ix)-this->n0/2.;
             const double ry=double(iy)-this->n1/2.;
             const double r = sqrt(pow(rx,2)+pow(ry,2));
             if(r>this->n0/4.){
                 if (point_up){
                     m(ix,iy,af::span,2) = -1.;
                 }
                 else {
                     m(ix,iy,af::span,2) = 1.;
                 }
            }
         }
     }
     return m;
}
