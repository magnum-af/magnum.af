#include "micro_nonequi_demag.hpp"

//Energy calculation
//Edemag=-mu0/2 integral(M . Hdemag) dx
double NonEquiDemagField::E(const State& state){
  return -constants::mu0/2. * material.ms * afvalue(sum(sum(sum(sum(h(state)*state.m,0),1),2),3)) * mesh.dx * mesh.dy * mesh.dz;
}

double NonEquiDemagField::E(const State& state, const af::array& h){
  return -constants::mu0/2. * material.ms * afvalue(sum(sum(sum(sum(h * state.m,0),1),2),3)) * mesh.dx * mesh.dy * mesh.dz;
}

void NonEquiDemagField::print_Nfft(){
    af::print("Nfft=", Nfft);
}

af::array N_cpp_alloc(int n0_exp, int n1_exp, int n2_exp, double dx, double dy, double dz);

NonEquiDemagField::NonEquiDemagField (Mesh meshin, Material paramin, bool verbose, bool caching, unsigned nthreads) : material(paramin),mesh(meshin), nthreads(nthreads > 0 ? nthreads : std::thread::hardware_concurrency()){
    af::timer demagtimer = af::timer::start();
    if (caching == false){
        Nfft=N_cpp_alloc(mesh.n0_exp,mesh.n1_exp,mesh.n2_exp,mesh.dx,mesh.dy,mesh.dz);
        if (verbose) printf("%s Starting Demag Tensor Assembly on %u out of %u threads.\n", Info(), this->nthreads, std::thread::hardware_concurrency());
        if (verbose) printf("time demag init [af-s]: %f\n", af::timer::stop(demagtimer));
    }
    else{
        std::string magafdir = setup_magafdir();
        unsigned long long maxsize = 2000000;
        unsigned long long reducedsize = 1000000;
        std::string nfft_id = "n0exp_"+std::to_string(mesh.n0_exp)+"_n1exp_"+std::to_string(mesh.n1_exp)+"_n2exp_"+std::to_string(mesh.n2_exp)+"_dx_"+std::to_string(1e9*mesh.dx)+"_dy_"+std::to_string(1e9*mesh.dy)+"_dz_"+std::to_string(1e9*mesh.dz);
        std::string path_to_nfft_cached = magafdir+nfft_id;
        int checkarray=-1;
        if (exists(path_to_nfft_cached)){
            try{
                checkarray = af::readArrayCheck(path_to_nfft_cached.c_str(), "");
            }
            catch (const af::exception& e){
                printf("Warning, af::readArrayCheck failed, omit reading array.\n%s\n", e.what());
            }
        }
        if(checkarray > -1){
            if (verbose) printf("Reading demag tensor from '%s'\n", path_to_nfft_cached.c_str());
            Nfft = af::readArray(path_to_nfft_cached.c_str(), "");
        }
        else{
            if (verbose) printf("%s Starting Demag Tensor Assembly on %u out of %u threads.\n", Info(), this->nthreads, std::thread::hardware_concurrency());
            Nfft=N_cpp_alloc(mesh.n0_exp,mesh.n1_exp,mesh.n2_exp,mesh.dx,mesh.dy,mesh.dz);
            unsigned long long magafdir_size_in_bytes = GetDirSize(magafdir);
            if (magafdir_size_in_bytes > maxsize){
                if (verbose) printf("Maintainance: size of '%s' is %f GB > %f GB, removing oldest files until size < %f GB\n", magafdir.c_str(), (double)magafdir_size_in_bytes/1e6, (double)maxsize/1e6, (double)reducedsize/1e6);
                remove_oldest_files_until_size(magafdir.c_str(), reducedsize, verbose);
                if (verbose) printf("Maintainance finished: '%s' has now %f GB\n", magafdir.c_str(), (double)GetDirSize(magafdir)/1e6);
            }
            if (GetDirSize(magafdir) < maxsize){
                try{
                    if (verbose) printf("Saving demag tensor to'%s'\n", path_to_nfft_cached.c_str());
                    af::saveArray("", Nfft, path_to_nfft_cached.c_str());
                    if (verbose) printf("Saved demag tensor to'%s'\n", path_to_nfft_cached.c_str());
                }
                catch (const af::exception& e){
                    printf("Warning, af::saveArray failed, omit saving demag tensor.\n%s\n", e.what());
                }
            }
        }
        if (verbose) printf("time demag init [af-s]: %f\n", af::timer::stop(demagtimer));
    }
}


af::array NonEquiDemagField::h(const State&  state){
    timer_demagsolve = af::timer::start();
  // FFT with zero-padding of the m field
  af::array mfft;
  if (mesh.n2_exp == 1){
      if (state.Ms.isempty()) mfft=af::fftR2C<2>(material.ms * state.m,af::dim4(mesh.n0_exp,mesh.n1_exp));
      else mfft=af::fftR2C<2>(state.Ms * state.m,af::dim4(mesh.n0_exp,mesh.n1_exp));
  }
  else {
      if (state.Ms.isempty()) mfft=af::fftR2C<3>(material.ms * state.m,af::dim4(mesh.n0_exp,mesh.n1_exp,mesh.n2_exp));
      else  mfft=af::fftR2C<3>(state.Ms * state.m,af::dim4(mesh.n0_exp,mesh.n1_exp,mesh.n2_exp));
  }

  // Pointwise product
  af::array hfft=af::array (mesh.n0_exp/2+1,mesh.n1_exp,mesh.n2_exp,3,c64);
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
  af::array h_field;
  if (mesh.n2_exp == 1){
    h_field=af::fftC2R<2>(hfft);
    if(material.afsync) af::sync();
    cpu_time += af::timer::stop(timer_demagsolve);
    return h_field(af::seq(0,mesh.n0_exp/2-1),af::seq(0,mesh.n1_exp/2-1));
  }
  else {
    h_field=af::fftC2R<3>(hfft);
    if(material.afsync) af::sync();
    cpu_time += af::timer::stop(timer_demagsolve);
    return h_field(af::seq(0,mesh.n0_exp/2-1),af::seq(0,mesh.n1_exp/2-1),af::seq(0,mesh.n2_exp/2-1),af::span);
  }
}

namespace newell{

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

    double F2(const double x, const double y, const double z){
        return newellf(x, y, z) - newellf(x, 0, z) - newellf(x, y, 0) + newellf(x, 0, 0);
    }

    double F1(const double x, const double y, const double z, const double dz, const double dZ){
        return F2(x, y, z + dZ) - F2(x, y, z) - F2(x, y, z - dz + dZ) + F2(x, y, z - dz);//TODO check dz vs dZ in first and last term
    }

    double F0(const double x, const double y, const double z, const double dy, const double dY, const double dz, const double dZ){
        return F1(x, y + dY, z, dz, dZ) - F1(x, y, z, dz, dZ) - F1(x, y - dy + dY, z, dz, dZ) + F1(x, y - dy, z, dz, dZ);//TODO check dz vs dZ in first and last term
    }
    double Nxx_nonequi(const int ix, const int iy, const int iz, const double dx, const double dy, const double dz, const double dX, const double dY, const double dZ){
        //TODO check def of xyz and tau
        const double x = dx*ix;
        const double y = dy*iy;
        const double z = dz*iz;
        const double tau = dx * dy * dz;//TODO check
        return -1./(4.0 * M_PI * tau) * ( \
                  F0(x          , y, z, dy, dY, dz, dZ) \
                - F0(x - dx     , y, z, dy, dY, dz, dZ) \
                - F0(x + dX     , y, z, dy, dY, dz, dZ) \
                + F0(x - dx + dX, y, z, dy, dY, dz, dZ));
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

}


struct LoopInfo {
    LoopInfo(){}
    LoopInfo(int i0_start, int i0_end, int n0_exp, int n1_exp, int n2_exp,  double dx,  double dy,  double dz):
        i0_start(i0_start), i0_end(i0_end), n0_exp(n0_exp), n1_exp(n1_exp), n2_exp(n2_exp), dx(dx), dy(dy), dz(dz){}
    int i0_start;
    int i0_end;
    int n0_exp;
    int n1_exp;
    int n2_exp;
    double dx;
    double dy;
    double dz;
};


double* N_nonequi_setup = NULL;


void* nonequi_setup_N(void* arg)
{
    LoopInfo* loopinfo = static_cast<LoopInfo*>(arg);
    for (int i0 = loopinfo->i0_start; i0 < loopinfo->i0_end; i0++){
        const int j0 = (i0 + loopinfo->n0_exp/2) % loopinfo->n0_exp - loopinfo->n0_exp/2;
        for (int i1 = 0; i1 < loopinfo->n1_exp; i1++){
            const int j1 = (i1 + loopinfo->n1_exp/2) % loopinfo->n1_exp - loopinfo->n1_exp/2;
            //for (int i2 = 0; i2 < loopinfo->n2_exp; i2++ ){
                //const int j2 = (i2 + loopinfo->n2_exp/2) % loopinfo->n2_exp - loopinfo->n2_exp/2;
                //const int idx = 6*(i2+loopinfo->n2_exp*(i1+loopinfo->n1_exp*i0));
                //std::cout << "j2=" << j2 << std::endl;
                const int idx = 6*(loopinfo->n2_exp*(i1+loopinfo->n1_exp*i0));
                const int j2 = 1;//TODO data structure for all j
                //TODO test with dx = dX, dy = dY, dz = dZ
                N_nonequi_setup[idx+0] = newell::Nxx_nonequi(j0, j1, j2, loopinfo->dx, loopinfo->dy, loopinfo->dz, loopinfo->dx, loopinfo->dy, loopinfo->dz);
                N_nonequi_setup[idx+1] = newell::Nxxg(j0, j1, j2, loopinfo->dx, loopinfo->dy, loopinfo->dz);
                N_nonequi_setup[idx+2] = newell::Nxxg(j0, j2, j1, loopinfo->dx, loopinfo->dz, loopinfo->dy);
                N_nonequi_setup[idx+3] = newell::Nxx_nonequi(j1, j2, j0, loopinfo->dy, loopinfo->dz, loopinfo->dx, loopinfo->dy, loopinfo->dz, loopinfo->dx);
                N_nonequi_setup[idx+4] = newell::Nxxg(j1, j2, j0, loopinfo->dy, loopinfo->dz, loopinfo->dx);
                N_nonequi_setup[idx+5] = newell::Nxx_nonequi(j2, j0, j1, loopinfo->dz, loopinfo->dx, loopinfo->dy, loopinfo->dz, loopinfo->dx, loopinfo->dy);

                //N_nonequi_setup[idx+0] = newell::Nxxf(j0, j1, j2, loopinfo->dx, loopinfo->dy, loopinfo->dz);
                //N_nonequi_setup[idx+1] = newell::Nxxg(j0, j1, j2, loopinfo->dx, loopinfo->dy, loopinfo->dz);
                //N_nonequi_setup[idx+2] = newell::Nxxg(j0, j2, j1, loopinfo->dx, loopinfo->dz, loopinfo->dy);
                //N_nonequi_setup[idx+3] = newell::Nxxf(j1, j2, j0, loopinfo->dy, loopinfo->dz, loopinfo->dx);
                //N_nonequi_setup[idx+4] = newell::Nxxg(j1, j2, j0, loopinfo->dy, loopinfo->dz, loopinfo->dx);
                //N_nonequi_setup[idx+5] = newell::Nxxf(j2, j0, j1, loopinfo->dz, loopinfo->dx, loopinfo->dy);
            //}
        }
    }
    return NULL;
} 


af::array NonEquiDemagField::N_cpp_alloc(int n0_exp, int n1_exp, int n2_exp, double dx, double dy, double dz){
    std::thread t[nthreads];
    struct LoopInfo loopinfo[nthreads];
    for (unsigned i = 0; i < nthreads; i++){
        unsigned start = i * (double)n0_exp/nthreads;
        unsigned end = (i +1) * (double)n0_exp/nthreads;
        loopinfo[i]=LoopInfo(start, end, n0_exp, n1_exp, n2_exp, dx, dy, dz);
    }

    N_nonequi_setup = new double[n0_exp*n1_exp*n2_exp*6];

    for (unsigned i = 0; i < nthreads; i++){
        t[i] = std::thread(nonequi_setup_N, &loopinfo[i]);
     }

    for (unsigned i = 0; i < nthreads; i++){
        t[i].join();
     }
    af::array Naf(6,n2_exp,n1_exp,n0_exp,N_nonequi_setup);
    Naf=af::reorder(Naf,3,2,1,0);
    delete [] N_nonequi_setup;
    N_nonequi_setup = NULL;
    af::print("Nonequi Naf", Naf);
    Naf = af::fftR2C<2>(Naf);
    af::print("Nonequi Nfft", Naf);
    //if (n2_exp == 1){
    //    Naf = af::fftR2C<2>(Naf);
    //}
    //else {
    //    Naf = af::fftR2C<3>(Naf);
    //}
    return Naf;
}
