#include "micro_demag_nonequi.hpp"

//Energy calculation
//Edemag=-mu0/2 integral(M . Hdemag) dx
double NonEquiDemagField::E(const State& state){
    return -constants::mu0/2. * state.material.ms * afvalue(sum(sum(sum(sum(h(state)*state.m,0),1),2),3)) * state.mesh.dx * state.mesh.dy * state.mesh.dz;
}

double NonEquiDemagField::E(const State& state, const af::array& h){
    return -constants::mu0/2. * state.material.ms * afvalue(sum(sum(sum(sum(h * state.m,0),1),2),3)) * state.mesh.dx * state.mesh.dy * state.mesh.dz;
}

void NonEquiDemagField::print_Nfft(){
    af::print("Nfft=", Nfft);
}


NonEquiDemagField::NonEquiDemagField (Mesh mesh, std::vector<double> z_spacing, bool verbose, bool caching, unsigned nthreads) : nthreads(nthreads > 0 ? nthreads : std::thread::hardware_concurrency()){
    af::timer demagtimer = af::timer::start();
    if (caching == false){
        Nfft=N_cpp_alloc(mesh.n0_exp, mesh.n1_exp, mesh.n2, mesh.dx, mesh.dy, mesh.dz, z_spacing);
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
            Nfft=N_cpp_alloc(mesh.n0_exp,mesh.n1_exp,mesh.n2_exp,mesh.dx,mesh.dy,mesh.dz, z_spacing);
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
    if (state.mesh.n2_exp == 1){
        if (state.Ms.isempty()) mfft=af::fftR2C<2>(state.material.ms * state.m,af::dim4(state.mesh.n0_exp,state.mesh.n1_exp));
        else mfft=af::fftR2C<2>(state.Ms * state.m,af::dim4(state.mesh.n0_exp,state.mesh.n1_exp));
    }
    else {
        if (state.Ms.isempty()) mfft=af::fftR2C<3>(state.material.ms * state.m,af::dim4(state.mesh.n0_exp,state.mesh.n1_exp,state.mesh.n2_exp));
        else  mfft=af::fftR2C<3>(state.Ms * state.m,af::dim4(state.mesh.n0_exp,state.mesh.n1_exp,state.mesh.n2_exp));
    }

    // Pointwise product
    af::array hfft=af::array (state.mesh.n0_exp/2+1, state.mesh.n1_exp, state.mesh.n2_exp, 3, c64);
    hfft(af::span, af::span, af::span, 0) = Nfft(af::span, af::span, af::span, 0) * mfft(af::span, af::span, af::span, 0)
                                          + Nfft(af::span, af::span, af::span, 1) * mfft(af::span, af::span, af::span, 1)
                                          + Nfft(af::span, af::span, af::span, 2) * mfft(af::span, af::span, af::span, 2);
    hfft(af::span, af::span, af::span, 1) = Nfft(af::span, af::span, af::span, 1) * mfft(af::span, af::span, af::span, 0)
                                          + Nfft(af::span, af::span, af::span, 3) * mfft(af::span, af::span, af::span, 1)
                                          + Nfft(af::span, af::span, af::span, 4) * mfft(af::span, af::span, af::span, 2);
    hfft(af::span, af::span, af::span, 2) = Nfft(af::span, af::span, af::span, 2) * mfft(af::span, af::span, af::span, 0)
                                          + Nfft(af::span, af::span, af::span, 4) * mfft(af::span, af::span, af::span, 1)
                                          + Nfft(af::span, af::span, af::span, 5) * mfft(af::span, af::span, af::span, 2);

    // IFFT reversing padding
    af::array h_field;
    if (state.mesh.n2_exp == 1){
        h_field=af::fftC2R<2>(hfft);
        if(state.material.afsync) af::sync();
        cpu_time += af::timer::stop(timer_demagsolve);
        return h_field(af::seq(0,state.mesh.n0_exp/2-1),af::seq(0,state.mesh.n1_exp/2-1));
    }
    else {
        h_field=af::fftC2R<3>(hfft);
        if(state.material.afsync) af::sync();
        cpu_time += af::timer::stop(timer_demagsolve);
        return h_field(af::seq(0,state.mesh.n0_exp/2-1),af::seq(0,state.mesh.n1_exp/2-1),af::seq(0,state.mesh.n2_exp/2-1),af::span);
    }
}

namespace newell_nonequi{

    double f(double x, double y, double z){
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
    
    double g(double x, double y, double z){
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

    //double F2(const double x, const double y, const double z){
    //    return f(x, y, z);
    //    //Last three terms cancel out//return f(x, y, z) - f(x, 0, z) - f(x, y, 0) + f(x, 0, 0);
    //}

    double F1(const double x, const double y, const double z, const double dz, const double dZ){
        return f(x, y, z + dZ) - f(x, y, z) - f(x, y, z - dz + dZ) + f(x, y, z - dz);//TODO check dz vs dZ in first and last term
    }

    double F0(const double x, const double y, const double z, const double dy, const double dY, const double dz, const double dZ){
        return F1(x, y + dY, z, dz, dZ) - F1(x, y, z, dz, dZ) - F1(x, y - dy + dY, z, dz, dZ) + F1(x, y - dy, z, dz, dZ);//TODO check dz vs dZ in first and last term
    }


    double Nxx_nonequi(const double x, const double y, const double z, const double dx, const double dy, const double dz, const double dX, const double dY, const double dZ){
        //TODO check def of xyz and tau
        const double tau = dx * dy * dz;//TODO check
        return -1./(4.0 * M_PI * tau) * ( \
                  F0(x          , y, z, dy, dY, dz, dZ) \
                - F0(x - dx     , y, z, dy, dY, dz, dZ) \
                - F0(x + dX     , y, z, dy, dY, dz, dZ) \
                + F0(x - dx + dX, y, z, dy, dY, dz, dZ));
    }

    //double G2(const double x, const double y, const double z){
    //    return g(x, y, z);
    //    //return g(x, y, z) - g(x, y, 0);
    //    //return g(x, y, z) - g(x, 0, z) - g(x, y, 0) + g(x, 0, 0);
    //}

    double G1(const double x, const double y, const double z, const double dz, const double dZ){
        return g(x, y, z + dZ) - g(x, y, z) - g(x, y, z - dz + dZ) + g(x, y, z - dz);//TODO check dz vs dZ in first and last term
    }

    double G0(const double x, const double y, const double z, const double dy, const double dY, const double dz, const double dZ){
        return G1(x, y + dY, z, dz, dZ) - G1(x, y, z, dz, dZ) - G1(x, y - dy + dY, z, dz, dZ) + G1(x, y - dy, z, dz, dZ);//TODO check dz vs dZ in first and last term
    }

    double Nxy_nonequi(const double x, const double y, const double z, const double dx, const double dy, const double dz, const double dX, const double dY, const double dZ){
        //TODO check def of xyz and tau
        const double tau = dx * dy * dz;//TODO check
        return -1./(4.0 * M_PI * tau) * ( \
                  G0(x          , y, z, dy, dY, dz, dZ) \
                - G0(x - dx     , y, z, dy, dY, dz, dZ) \
                - G0(x + dX     , y, z, dy, dY, dz, dZ) \
                + G0(x - dx + dX, y, z, dy, dY, dz, dZ));
    }

    class NonequiLoopInfo {
        public:
        NonequiLoopInfo(){}
        NonequiLoopInfo(int i0_start, int i0_end, int n0_exp, int n1_exp, int n2,  double dx,  double dy,  double dz):
            i0_start(i0_start), i0_end(i0_end), n0_exp(n0_exp), n1_exp(n1_exp), n2(n2), dx(dx), dy(dy), dz(dz){}
        int i0_start;
        int i0_end;
        int n0_exp;
        int n1_exp;
        int n2;
        double dx;
        double dy;
        double dz;
        static std::vector<double> z_spacing;
    };

    std::vector<double> newell_nonequi::NonequiLoopInfo::z_spacing; //Initialize static member

    double nonequi_index_distance(const std::vector<double> spacings, const unsigned i, const unsigned j){
        //Calculates the distance beween indices i and j
        //TODO check alternative definition, define way to count
        //Now: 1_|2__|3_
        //Now: °_|°__|°_
        //1:2   _
        //1:3   _  __
        //2:3      __
        //but dz3 is never used, TODO check!
        // Examples describing definition:
        // d(i,i) = 0
        // d(i,i+1) =  z_spacing[i]
        // d(i,i-1) = -z_spacing[i]
        // d(i,i+2) =  z_spacing[i] + z_spacing[i+1]
    
        double result = 0;
        if(i > j){
            for (unsigned k = i; k > j; k--){
                result -= spacings.at(k-1);//TODO taking minus_index value
                //std::cout << "k-1=" << k-1 << "  " << spacings[k] << "  " << result <<  std::endl;
            }
        }
        else {
            for (unsigned k = i; k < j; k++){
                result += spacings.at(k);
                //std::cout << k << "  " << spacings[k] << "  " << result <<  std::endl;
            }
        }
        return result;
    }

    
    double* N_nonequi_setup = NULL;
    
    void* nonequi_setup_N(void* arg)
    {
        newell_nonequi::NonequiLoopInfo* loopinfo = static_cast<newell_nonequi::NonequiLoopInfo*>(arg);
        for (int i0 = loopinfo->i0_start; i0 < loopinfo->i0_end; i0++){
            const int j0 = (i0 + loopinfo->n0_exp/2) % loopinfo->n0_exp - loopinfo->n0_exp/2;
            for (int i1 = 0; i1 < loopinfo->n1_exp; i1++){
                const int j1 = (i1 + loopinfo->n1_exp/2) % loopinfo->n1_exp - loopinfo->n1_exp/2;
                for (int i2 = 0; i2 < loopinfo->n2; i2++ ){
                    for (int i3 = 0; i3 < loopinfo->n2; i3++ ){
                        //const int j2 = (i2 + loopinfo->n2/2) % loopinfo->n2 - loopinfo->n2/2;
                        //const int idx = 6*(i2+loopinfo->n2*(i1+loopinfo->n1_exp*i0));
                        const int idx = 6*(i3+loopinfo->n2*(i2+loopinfo->n2*(i1+loopinfo->n1_exp*i0)));
    
                        const double x = loopinfo->dx * (double)j0;
                        const double y = loopinfo->dy * (double)j1;
                        const double z = nonequi_index_distance(loopinfo->z_spacing, i2, i3);//loopinfo->dz * (double)j2;
                        //const int j2 = 1;//TODO data structure for all j
                        //std::cout << "j2=" << j2 << std::endl;
                        //const int idx = 6*(loopinfo->n2*(i1+loopinfo->n1_exp*i0));
                        //TODO test with dx = dX, dy = dY, dz = dZ
    
                        //TODO: check i2 == dz
                        //TODO: check i3 == dZ
                        //std::cout << i0 <<  " " << i1 <<  " " << i2 <<  " " << i3 <<  " " << z <<  " " << std::endl;
                        //std::cout << x << " "  << y << " " << z << std::endl;
                        //std::cout << z << std::endl;
                        newell_nonequi::N_nonequi_setup[idx+0] = newell_nonequi::Nxx_nonequi(x, y, z, loopinfo->dx, loopinfo->dy, loopinfo->z_spacing[i2], loopinfo->dx, loopinfo->dy, loopinfo->z_spacing[i3]);
                        newell_nonequi::N_nonequi_setup[idx+1] = newell_nonequi::Nxy_nonequi(x, y, z, loopinfo->dx, loopinfo->dy, loopinfo->z_spacing[i2], loopinfo->dx, loopinfo->dy, loopinfo->z_spacing[i3]);
                        newell_nonequi::N_nonequi_setup[idx+2] = newell_nonequi::Nxy_nonequi(x, z, y, loopinfo->dx, loopinfo->z_spacing[i2], loopinfo->dy, loopinfo->dx, loopinfo->z_spacing[i3], loopinfo->dy);
                        newell_nonequi::N_nonequi_setup[idx+3] = newell_nonequi::Nxx_nonequi(y, z, x, loopinfo->dy, loopinfo->z_spacing[i2], loopinfo->dx, loopinfo->dy, loopinfo->z_spacing[i3], loopinfo->dx);
                        newell_nonequi::N_nonequi_setup[idx+4] = newell_nonequi::Nxy_nonequi(y, z, x, loopinfo->dy, loopinfo->z_spacing[i2], loopinfo->dx, loopinfo->dy, loopinfo->z_spacing[i3], loopinfo->dx);
                        newell_nonequi::N_nonequi_setup[idx+5] = newell_nonequi::Nxx_nonequi(z, x, y, loopinfo->z_spacing[i2], loopinfo->dx, loopinfo->dy, loopinfo->z_spacing[i3], loopinfo->dx, loopinfo->dy);
                    }
                }
            }
        }
        return NULL;
    } 
}


af::array NonEquiDemagField::N_cpp_alloc(int n0_exp, int n1_exp, int n2, double dx, double dy, double dz, const std::vector<double> z_spacing){
    std::thread t[nthreads];
    struct newell_nonequi::NonequiLoopInfo loopinfo[nthreads];
    newell_nonequi::NonequiLoopInfo::z_spacing = z_spacing;
    for (unsigned i = 0; i < nthreads; i++){
        unsigned start = i * (double)n0_exp/nthreads;
        unsigned end = (i +1) * (double)n0_exp/nthreads;
        loopinfo[i]=newell_nonequi::NonequiLoopInfo(start, end, n0_exp, n1_exp, n2, dx, dy, dz);
    }

    newell_nonequi::N_nonequi_setup = new double[n0_exp * n1_exp * n2 * n2 * 6];

    for (unsigned i = 0; i < nthreads; i++){
        t[i] = std::thread(newell_nonequi::nonequi_setup_N, &loopinfo[i]);
     }

    for (unsigned i = 0; i < nthreads; i++){
        t[i].join();
     }
    af::array Naf(6, n2 * n2, n1_exp, n0_exp, newell_nonequi::N_nonequi_setup);
    Naf=af::reorder(Naf,3,2,1,0);
    todel_N = Naf;//TODO todel
    delete [] newell_nonequi::N_nonequi_setup;
    newell_nonequi::N_nonequi_setup = NULL;
    //af::print("Nonequi Naf", Naf);
    //TODO//Naf = af::fftR2C<2>(Naf);
    //af::print("Nonequi Nfft", Naf);

    //if (n2_exp == 1){
    std::cout << "Naf.dims() = " << Naf.dims() << std::endl;
    Naf = af::fftR2C<2>(Naf);
    //}
    //else {
    //    Naf = af::fftR2C<3>(Naf);
    //}
    return Naf;
}
