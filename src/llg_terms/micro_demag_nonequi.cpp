#include "micro_demag_nonequi.hpp"
#include "../misc.hpp"
#include "../func.hpp"
#include <iostream>
#include <vector>
#include <thread>

namespace magnumafcpp{


//Energy calculation: Edemag = - mu0/2 * integral(M . Hdemag) dx
double NonEquiDemagField::E(const State& state){
    return - constants::mu0/2. * state.integral_nonequimesh(h(state) * state.m);
}


double NonEquiDemagField::E(const State& state, const af::array& h){
    return - constants::mu0/2. * state.integral_nonequimesh(h * state.m);
}


NonEquiDemagField::NonEquiDemagField (NonequispacedMesh mesh, bool verbose, bool caching, unsigned nthreads) : nthreads(nthreads > 0 ? nthreads : std::thread::hardware_concurrency()){
    af::timer demagtimer = af::timer::start();
    if (caching == false){
        Nfft=calculate_N(mesh.nx_expanded, mesh.ny_expanded, mesh.nz, mesh.dx, mesh.dy, mesh.z_spacing);
        if (verbose) printf("%s Starting Demag Tensor Assembly on %u out of %u threads.\n", Info(), this->nthreads, std::thread::hardware_concurrency());
        if (verbose) printf("time demag init [af-s]: %f\n", af::timer::stop(demagtimer));
    }
    else{
        std::string magafdir = setup_magafdir();
        unsigned long long maxsize = 2000000;
        unsigned long long reducedsize = 1000000;
        //+std::to_string(1e9*mesh.dz);
        std::string dz_string;
        for (auto const& dz : mesh.z_spacing){
            dz_string.append(std::to_string(1e9 * dz));
        }
        std::string nfft_id = "n0exp_"+std::to_string(mesh.nx_expanded)+"_n1exp_"+std::to_string(mesh.ny_expanded)+"_nz_"+std::to_string(mesh.nz)+"_dx_"+std::to_string(1e9*mesh.dx)+"_dy_"+std::to_string(1e9*mesh.dy)+"_dz_"+dz_string;

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
            Nfft=calculate_N(mesh.nx_expanded, mesh.ny_expanded, mesh.nz, mesh.dx, mesh.dy, mesh.z_spacing);
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
    if (state.Ms_field.isempty()){
        mfft = af::fftR2C<2>(state.Ms * state.m, af::dim4(state.nonequimesh.nx_expanded, state.nonequimesh.ny_expanded));
    }
    else {
        mfft = af::fftR2C<2>(state.Ms_field * state.m, af::dim4(state.nonequimesh.nx_expanded, state.nonequimesh.ny_expanded));
    }

    // Pointwise product
    af::array hfft = af::constant(0.0, state.nonequimesh.nx_expanded/2+1, state.nonequimesh.ny_expanded, state.nonequimesh.nz, 3, c64);
    for (int i_source = 0; i_source < state.nonequimesh.nz; i_source++ ){
        for (int i_target = 0; i_target < state.nonequimesh.nz; i_target++ ){

            int zindex = util::ij2k(i_source, i_target, state.nonequimesh.nz);
            af::array nfft;
            if (i_source <= i_target){// This reflects the data structure of newell_nonequi::N_ptr. "<=" choosen by definition, could also be ">(=)" or just "<" when reconsidering N_ptr
                nfft = Nfft(af::span, af::span, zindex, af::span);
            }
            else {
                //swap indices to acces symmetric element ij -> ji
                nfft = Nfft(af::span, af::span, zindex, af::span);
                nfft (af::span, af::span, af::span, 2) = -1 * Nfft(af::span, af::span, zindex, 2);
                nfft (af::span, af::span, af::span, 4) = -1 * Nfft(af::span, af::span, zindex, 4);
            }

            hfft(af::span, af::span, i_target, 0) += nfft(af::span, af::span, af::span, 0) * mfft(af::span, af::span, i_source, 0)
                                                   + nfft(af::span, af::span, af::span, 1) * mfft(af::span, af::span, i_source, 1)
                                                   + nfft(af::span, af::span, af::span, 2) * mfft(af::span, af::span, i_source, 2);
            hfft(af::span, af::span, i_target, 1) += nfft(af::span, af::span, af::span, 1) * mfft(af::span, af::span, i_source, 0)
                                                   + nfft(af::span, af::span, af::span, 3) * mfft(af::span, af::span, i_source, 1)
                                                   + nfft(af::span, af::span, af::span, 4) * mfft(af::span, af::span, i_source, 2);
            hfft(af::span, af::span, i_target, 2) += nfft(af::span, af::span, af::span, 2) * mfft(af::span, af::span, i_source, 0)
                                                   + nfft(af::span, af::span, af::span, 4) * mfft(af::span, af::span, i_source, 1)
                                                   + nfft(af::span, af::span, af::span, 5) * mfft(af::span, af::span, i_source, 2);
        }
    }

    af::array one_over_tau_vec = af::array(1, 1, state.nonequimesh.nz, 1, f64);
    for (int i = 0; i < state.nonequimesh.nz; i++){
        one_over_tau_vec(0, 0, i, 0) = 1./(state.nonequimesh.dx * state.nonequimesh.dy * state.nonequimesh.z_spacing[i]);
        //std::cout << afvalue(one_over_tau_vec(0, 0, i, 0)) << "\n";
    }
    //af::print("tau", one_over_tau_vec);
    one_over_tau_vec = af::tile(one_over_tau_vec, state.nonequimesh.nx, state.nonequimesh.ny, 1, 3);
    //af::print("tau", one_over_tau_vec);

    // IFFT reversing padding
    af::array h_field;
    h_field=af::fftC2R<2>(hfft);
    if(state.material.afsync) af::sync();
    cpu_time += af::timer::stop(timer_demagsolve);
    //return h_field(af::seq(0, state.nonequimesh.nx_expanded/2-1), af::seq(0, state.nonequimesh.ny_expanded/2-1));
    return one_over_tau_vec * h_field(af::seq(0, state.nonequimesh.nx_expanded/2-1), af::seq(0, state.nonequimesh.ny_expanded/2-1));
}

namespace newell_nonequi{

    double f(double x, double y, double z){
      x=fabs(x);
      y=fabs(y);
      z=fabs(z);
      const double R = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
      const double xx = pow(x, 2);
      const double yy = pow(y, 2);
      const double zz = pow(z, 2);

      double result = 1.0 / 6.0 * (2.0*xx - yy - zz) * R;
      if(xx + zz > 0) result += y / 2.0 * (zz - xx) * asinh(y / (sqrt(xx + zz)));
      if(xx + yy > 0) result += z / 2.0 * (yy - xx) * asinh(z / (sqrt(xx + yy)));
      if(x  *  R > 0) result += - x*y*z * atan(y*z / (x * R));
      return result;
    }

    double g(double x, double y, double z){
      z=fabs(z);
      const double R = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
      const double xx = pow(x, 2);
      const double yy = pow(y, 2);
      const double zz = pow(z, 2);

      double result = - x*y * R / 3.0;
      if(xx + yy > 0) result += x*y*z * asinh(z / (sqrt(xx + yy)));
      if(yy + zz > 0) result += y / 6.0 * (3.0 * zz - yy) * asinh(x / (sqrt(yy + zz)));
      if(xx + zz > 0) result += x / 6.0 * (3.0 * zz - xx) * asinh(y / (sqrt(xx + zz)));
      if(z  *  R > 0) result += - pow(z, 3) / 6.0     * atan(x*y / (z * R));
      if(y  *  R!= 0) result += - z * yy / 2.0 * atan(x*z / (y * R));
      if(x  *  R!= 0) result += - z * xx / 2.0 * atan(y*z / (x * R));
      return result;
    }

    //double F2(const double x, const double y, const double z){
    //    return f(x, y, z);
    //    //Last three terms cancel out//return f(x, y, z) - f(x, 0, z) - f(x, y, 0) + f(x, 0, 0);
    //}

    double F1(const double x, const double y, const double z, const double dz, const double dZ){
        return f(x, y, z + dZ) - f(x, y, z) - f(x, y, z - dz + dZ) + f(x, y, z - dz);
    }

    double F0(const double x, const double y, const double z, const double dy, const double dY, const double dz, const double dZ){
        return F1(x, y + dY, z, dz, dZ) - F1(x, y, z, dz, dZ) - F1(x, y - dy + dY, z, dz, dZ) + F1(x, y - dy, z, dz, dZ);
    }


    double Nxx(const double x, const double y, const double z, const double dx, const double dy, const double dz, const double dX, const double dY, const double dZ){
        // x, y, z is vector from source cuboid to target cuboid
        // dx, dy, dz are source cuboid dimensions
        // dX, dY, dZ are target cuboid dimensions
        //const double tau = dX * dY * dZ;// Defining dX, dY, dZ as target cuboid (one could alternatively choose dx, dy, dz with implications on x, y, z)
        //return -1./(4.0 * M_PI * tau) * (
        return -1./(4.0 * M_PI) * ( \
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
        return g(x, y, z + dZ) - g(x, y, z) - g(x, y, z - dz + dZ) + g(x, y, z - dz);
    }

    double G0(const double x, const double y, const double z, const double dy, const double dY, const double dz, const double dZ){
        return G1(x, y + dY, z, dz, dZ) - G1(x, y, z, dz, dZ) - G1(x, y - dy + dY, z, dz, dZ) + G1(x, y - dy, z, dz, dZ);
    }

    double Nxy(const double x, const double y, const double z, const double dx, const double dy, const double dz, const double dX, const double dY, const double dZ){
        // x, y, z is vector from source cuboid to target cuboid
        // dx, dy, dz are source cuboid dimensions
        // dX, dY, dZ are target cuboid dimensions
        //const double tau = dX * dY * dZ;// Defining dX, dY, dZ as target cuboid (one could alternatively choose dx, dy, dz with implications on x, y, z)
        //return -1./(4.0 * M_PI * tau) * (
        return -1./(4.0 * M_PI) * ( \
                  G0(x          , y, z, dy, dY, dz, dZ) \
                - G0(x - dx     , y, z, dy, dY, dz, dZ) \
                - G0(x + dX     , y, z, dy, dY, dz, dZ) \
                + G0(x - dx + dX, y, z, dy, dY, dz, dZ));
    }

    class NonequiLoopInfo {
        public:
        NonequiLoopInfo(){}
        NonequiLoopInfo(int ix_start, int ix_end, int n0_exp, int n1_exp, int n2,  double dx,  double dy):
            ix_start(ix_start), ix_end(ix_end), n0_exp(n0_exp), n1_exp(n1_exp), n2(n2), dx(dx), dy(dy){}
        int ix_start;
        int ix_end;
        int n0_exp;
        int n1_exp;
        int n2;
        double dx;
        double dy;
        static std::vector<double> z_spacing;
    };

    std::vector<double> newell_nonequi::NonequiLoopInfo::z_spacing; //Declare static member

    double nonequi_index_distance(const std::vector<double> spacings, const unsigned i, const unsigned j, const bool verbose){
        //Calculates the signed distance beween elements by summing up i < j: sum_(k=i)^(j-1)[spacings[k]] or i > j: sum_(k=j)^(i-1)[ - spacings[k]]
        //Note that spacings[spacings.size()-1] is never used
        if (verbose and (i == spacings.size() or j == spacings.size())){
            printf("%s in nonequi_index_distance: index == vector.size(), the distance includes the last element which is not wanted behaviour\n", Warning());
        }

        double result = 0;
        if(i > j){
            for (unsigned k = i; k > j; k--){
                result -= spacings.at(k-1);
            }
        }
        else {
            for (unsigned k = i; k < j; k++){
                result += spacings.at(k);
            }
        }
        return result;
    }


    double* N_ptr = NULL;

    void* init_N(void* arg)
    {
        newell_nonequi::NonequiLoopInfo* loopinfo = static_cast<newell_nonequi::NonequiLoopInfo*>(arg);
        for (int ix = loopinfo->ix_start; ix < loopinfo->ix_end; ix++){
            const int jx = (ix + loopinfo->n0_exp/2) % loopinfo->n0_exp - loopinfo->n0_exp/2;
            for (int iy = 0; iy < loopinfo->n1_exp; iy++){
                const int jy = (iy + loopinfo->n1_exp/2) % loopinfo->n1_exp - loopinfo->n1_exp/2;
                for (int i_source = 0; i_source < loopinfo->n2; i_source++ ){
                    for (int i_target = 0; i_target < loopinfo->n2; i_target++ ){

                        if (i_source <= i_target){
                            const int idx = 6*(util::ij2k(i_source, i_target, loopinfo->n2) + ((loopinfo->n2 * (loopinfo->n2 + 1))/2)*(iy+loopinfo->n1_exp*ix));
                            //std::cout << "idx=" << idx << " of " << loopinfo->n0_exp * loopinfo->n1_exp * (loopinfo->n2 * (loopinfo->n2 + 1))/2 * 6 << std::endl;
                            const double x = loopinfo->dx * (double)jx;
                            const double y = loopinfo->dy * (double)jy;
                            const double z = nonequi_index_distance(loopinfo->z_spacing, i_source, i_target);

                            newell_nonequi::N_ptr[idx+0] = newell_nonequi::Nxx(x, y, z, loopinfo->dx, loopinfo->dy, loopinfo->z_spacing[i_source], loopinfo->dx, loopinfo->dy, loopinfo->z_spacing[i_target]);
                            newell_nonequi::N_ptr[idx+1] = newell_nonequi::Nxy(x, y, z, loopinfo->dx, loopinfo->dy, loopinfo->z_spacing[i_source], loopinfo->dx, loopinfo->dy, loopinfo->z_spacing[i_target]);
                            newell_nonequi::N_ptr[idx+2] = newell_nonequi::Nxy(x, z, y, loopinfo->dx, loopinfo->z_spacing[i_source], loopinfo->dy, loopinfo->dx, loopinfo->z_spacing[i_target], loopinfo->dy);
                            newell_nonequi::N_ptr[idx+3] = newell_nonequi::Nxx(y, z, x, loopinfo->dy, loopinfo->z_spacing[i_source], loopinfo->dx, loopinfo->dy, loopinfo->z_spacing[i_target], loopinfo->dx);
                            newell_nonequi::N_ptr[idx+4] = newell_nonequi::Nxy(y, z, x, loopinfo->dy, loopinfo->z_spacing[i_source], loopinfo->dx, loopinfo->dy, loopinfo->z_spacing[i_target], loopinfo->dx);
                            newell_nonequi::N_ptr[idx+5] = newell_nonequi::Nxx(z, x, y, loopinfo->z_spacing[i_source], loopinfo->dx, loopinfo->dy, loopinfo->z_spacing[i_target], loopinfo->dx, loopinfo->dy);
                        }
                    }
                }
            }
        }
        return NULL;
    }
}


af::array NonEquiDemagField::calculate_N(int n0_exp, int n1_exp, int n2, double dx, double dy, const std::vector<double> z_spacing){
    std::vector<std::thread> t;
    std::vector<newell_nonequi::NonequiLoopInfo> loopinfo;
    newell_nonequi::NonequiLoopInfo::z_spacing = z_spacing;
    for (unsigned i = 0; i < nthreads; i++){
        unsigned start = i * (double)n0_exp/nthreads;
        unsigned end = (i +1) * (double)n0_exp/nthreads;
        loopinfo.push_back(newell_nonequi::NonequiLoopInfo(start, end, n0_exp, n1_exp, n2, dx, dy));
    }

    newell_nonequi::N_ptr = new double[n0_exp * n1_exp * (n2 * (n2 + 1))/2 * 6];

    for (unsigned i = 0; i < nthreads; i++){
        t.push_back( std::thread(newell_nonequi::init_N, &loopinfo[i]) );
     }

    for (unsigned i = 0; i < nthreads; i++){
        t[i].join();
     }
    af::array Naf(6, (n2 * (n2 + 1))/2, n1_exp, n0_exp, newell_nonequi::N_ptr);
    Naf=af::reorder(Naf, 3, 2, 1, 0);
    delete [] newell_nonequi::N_ptr;
    newell_nonequi::N_ptr = NULL;
    Naf = af::fftR2C<2>(Naf);
    return Naf;
}
}// namespace magnumafcpp
