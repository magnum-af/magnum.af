#include "demag.hpp"

//Energy calculation
//Edemag=-mu0/2 integral(M . Hdemag) dx
double DemagSolver::E(const State& state){
  return -param.mu0/2. * param.ms * afvalue(sum(sum(sum(sum(h(state)*state.m,0),1),2),3)) * mesh.dx * mesh.dy * mesh.dz; 
}

af::array N_cpp_alloc(int n0_exp, int n1_exp, int n2_exp, double dx, double dy, double dz);

DemagSolver::DemagSolver (Mesh meshin, Param paramin) : param(paramin),mesh(meshin){
  Nfft=N_cpp_alloc(mesh.n0_exp,mesh.n1_exp,mesh.n2_exp,mesh.dx,mesh.dy,mesh.dz);
  h_field   =af::array (mesh.n0_exp    ,mesh.n1_exp,mesh.n2_exp,3,f64);
  hfft=af::array (mesh.n0_exp/2+1,mesh.n1_exp,mesh.n2_exp,3,c64);
}


af::array DemagSolver::h(const State&  state){
    timer_demagsolve = af::timer::start();
  // FFT with zero-padding of the m field
  if (mesh.n2_exp == 1){
    mfft=af::fftR2C<2>(state.m,af::dim4(mesh.n0_exp,mesh.n1_exp));
  }
  else {
    mfft=af::fftR2C<3>(state.m,af::dim4(mesh.n0_exp,mesh.n1_exp,mesh.n2_exp));
  }

  // Pointwise product
  hfft(af::span,af::span,af::span,0)= Nfft(af::span,af::span,af::span,0) * mfft(af::span,af::span,af::span,0) 
                                 + Nfft(af::span,af::span,af::span,1) * mfft(af::span,af::span,af::span,1)
                                 + Nfft(af::span,af::span,af::span,2) * mfft(af::span,af::span,af::span,2);
  hfft(af::span,af::span,af::span,1)= Nfft(af::span,af::span,af::span,1) * mfft(af::span,af::span,af::span,0) 
                                 + Nfft(af::span,af::span,af::span,3) * mfft(af::span,af::span,af::span,1)
                                 + Nfft(af::span,af::span,af::span,4) * mfft(af::span,af::span,af::span,2);
  hfft(af::span,af::span,af::span,2)= Nfft(af::span,af::span,af::span,2) * mfft(af::span,af::span,af::span,0)
                                 + Nfft(af::span,af::span,af::span,4) * mfft(af::span,af::span,af::span,1)
                                 + Nfft(af::span,af::span,af::span,5) * mfft(af::span,af::span,af::span,2);

  // IFFT reversing padding
  if (mesh.n2_exp == 1){
    h_field=af::fftC2R<2>(hfft);
    if(param.afsync) af::sync();
    cpu_time += af::timer::stop(timer_demagsolve);
    return param.ms * h_field(af::seq(0,mesh.n0_exp/2-1),af::seq(0,mesh.n1_exp/2-1));
  }
  else {
    h_field=af::fftC2R<3>(hfft);
    if(param.afsync) af::sync();
    cpu_time += af::timer::stop(timer_demagsolve);
    return param.ms * h_field(af::seq(0,mesh.n0_exp/2-1),af::seq(0,mesh.n1_exp/2-1),af::seq(0,mesh.n2_exp/2-1),af::span);
  }
}

double newellf(double x, double y, double z){
  x=fabs(x);
  y=fabs(y);
  z=fabs(z);
  const double R = sqrt(pow(x,2) + pow(y,2) + pow(z,2));
  const double xx = pow(x,2);
  const double yy = pow(y,2);
  const double zz = pow(z,2);

  double result = 1.0 / 6.0 * (2.0*xx - yy - zz) * R;
  if(xx + zz > 0) result += y / 2.0 * (zz - xx) * asinh(y / (sqrt(xx + zz)));
  if(xx + yy > 0) result += z / 2.0 * (yy - xx) * asinh(z / (sqrt(xx + yy)));
  if(x  *  R > 0) result += - x*y*z * atan(y*z / (x * R));
  return result;
}

double newellg(double x, double y, double z){
  z=fabs(z);
  const double R = sqrt(pow(x,2) + pow(y,2) + pow(z,2));
  const double xx = pow(x,2);
  const double yy = pow(y,2);
  const double zz = pow(z,2);

  double result = - x*y * R / 3.0;
  if(xx + yy > 0) result += x*y*z * asinh(z / (sqrt(xx + yy)));
  if(yy + zz > 0) result += y / 6.0 * (3.0 * zz - yy) * asinh(x / (sqrt(yy + zz)));
  if(xx + zz > 0) result += x / 6.0 * (3.0 * zz - xx) * asinh(y / (sqrt(xx + zz)));
  if(z  *  R > 0) result += - pow(z,3) / 6.0     * atan(x*y / (z * R));
  if(y  *  R!= 0) result += - z * yy / 2.0 * atan(x*z / (y * R));
  if(x  *  R!= 0) result += - z * xx / 2.0 * atan(y*z / (x * R));
  return result;
}

double Nxxf(int ix, int iy, int iz, double dx, double dy, double dz){
  double x = dx*ix;
  double y = dy*iy;
  double z = dz*iz;
  double result = 8.0 * newellf( x,    y,    z   ) \
         - 4.0 * newellf( x+dx, y,    z   ) \
         - 4.0 * newellf( x-dx, y,    z   ) \
         - 4.0 * newellf( x,    y+dy, z   ) \
         - 4.0 * newellf( x,    y-dy, z   ) \
         - 4.0 * newellf( x,    y   , z+dz) \
         - 4.0 * newellf( x,    y   , z-dz) \
         + 2.0 * newellf( x+dx, y+dy, z   ) \
         + 2.0 * newellf( x+dx, y-dy, z   ) \
         + 2.0 * newellf( x-dx, y+dy, z   ) \
         + 2.0 * newellf( x-dx, y-dy, z   ) \
         + 2.0 * newellf( x+dx, y   , z+dz) \
         + 2.0 * newellf( x+dx, y   , z-dz) \
         + 2.0 * newellf( x-dx, y   , z+dz) \
         + 2.0 * newellf( x-dx, y   , z-dz) \
         + 2.0 * newellf( x   , y+dy, z+dz) \
         + 2.0 * newellf( x   , y+dy, z-dz) \
         + 2.0 * newellf( x   , y-dy, z+dz) \
         + 2.0 * newellf( x   , y-dy, z-dz) \
         - 1.0 * newellf( x+dx, y+dy, z+dz) \
         - 1.0 * newellf( x+dx, y+dy, z-dz) \
         - 1.0 * newellf( x+dx, y-dy, z+dz) \
         - 1.0 * newellf( x+dx, y-dy, z-dz) \
         - 1.0 * newellf( x-dx, y+dy, z+dz) \
         - 1.0 * newellf( x-dx, y+dy, z-dz) \
         - 1.0 * newellf( x-dx, y-dy, z+dz) \
         - 1.0 * newellf( x-dx, y-dy, z-dz);
  return - result / (4.0 * M_PI * dx * dy * dz);
}

double Nxxg(int ix, int iy, int iz, double dx, double dy, double dz){
  double x = dx*ix;
  double y = dy*iy;
  double z = dz*iz;
  double result = 8.0 * newellg( x,    y,    z   ) \
                - 4.0 * newellg( x+dx, y,    z   ) \
                - 4.0 * newellg( x-dx, y,    z   ) \
                - 4.0 * newellg( x,    y+dy, z   ) \
                - 4.0 * newellg( x,    y-dy, z   ) \
                - 4.0 * newellg( x,    y   , z+dz) \
                - 4.0 * newellg( x,    y   , z-dz) \
                + 2.0 * newellg( x+dx, y+dy, z   ) \
                + 2.0 * newellg( x+dx, y-dy, z   ) \
                + 2.0 * newellg( x-dx, y+dy, z   ) \
                + 2.0 * newellg( x-dx, y-dy, z   ) \
                + 2.0 * newellg( x+dx, y   , z+dz) \
                + 2.0 * newellg( x+dx, y   , z-dz) \
                + 2.0 * newellg( x-dx, y   , z+dz) \
                + 2.0 * newellg( x-dx, y   , z-dz) \
                + 2.0 * newellg( x   , y+dy, z+dz) \
                + 2.0 * newellg( x   , y+dy, z-dz) \
                + 2.0 * newellg( x   , y-dy, z+dz) \
                + 2.0 * newellg( x   , y-dy, z-dz) \
                - 1.0 * newellg( x+dx, y+dy, z+dz) \
                - 1.0 * newellg( x+dx, y+dy, z-dz) \
                - 1.0 * newellg( x+dx, y-dy, z+dz) \
                - 1.0 * newellg( x+dx, y-dy, z-dz) \
                - 1.0 * newellg( x-dx, y+dy, z+dz) \
                - 1.0 * newellg( x-dx, y+dy, z-dz) \
                - 1.0 * newellg( x-dx, y-dy, z+dz) \
                - 1.0 * newellg( x-dx, y-dy, z-dz);
  result = - result / (4.0 * M_PI * dx * dy * dz);
  return result;
}

af::array N_cpp_alloc(int n0_exp, int n1_exp, int n2_exp, double dx, double dy, double dz){
  double* N = NULL;
  N = new double[n0_exp*n1_exp*n2_exp*6];
  for (int i0 = 0; i0 < n0_exp; i0++){
    const int j0 = (i0 + n0_exp/2) % n0_exp - n0_exp/2;
    for (int i1 = 0; i1 < n1_exp; i1++){
      const int j1 = (i1 + n1_exp/2) % n1_exp - n1_exp/2;
      for (int i2 = 0; i2 < n2_exp; i2++ ){
        const int j2 = (i2 + n2_exp/2) % n2_exp - n2_exp/2;
        const int idx = 6*(i2+n2_exp*(i1+n1_exp*i0));
        N[idx+0] = Nxxf(j0, j1, j2, dx, dy, dz);
        N[idx+1] = Nxxg(j0, j1, j2, dx, dy, dz);
        N[idx+2] = Nxxg(j0, j2, j1, dx, dz, dy);
        N[idx+3] = Nxxf(j1, j2, j0, dy, dz, dx);
        N[idx+4] = Nxxg(j1, j2, j0, dy, dz, dx);
        N[idx+5] = Nxxf(j2, j0, j1, dz, dx, dy);
      }
    }
  }
  af::array Naf(6,n2_exp,n1_exp,n0_exp,N);
  Naf=af::reorder(Naf,3,2,1,0);
  delete [] N;
  N = NULL;
  if (n2_exp == 1){
    Naf = af::fftR2C<2>(Naf);
  }
  else {
    Naf = af::fftR2C<3>(Naf);
  }
  return Naf;
}
