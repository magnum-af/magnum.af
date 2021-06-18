#include "nonequi/nonequi_demag_field.hpp"
#include "util/cache_manager.hpp"
#include "util/util.hpp"
#include "util/color_string.hpp"
#include "util/prime_factors.hpp"
#include <iostream>
#include <thread>
#include <vector>

namespace magnumafcpp {

/// Expanded cell size for demag FFT
inline unsigned nx_expanded(unsigned nx) { return 2 * nx; }
inline unsigned ny_expanded(unsigned ny) { return 2 * ny; }

af::array NonequiDemagField::impl_H_in_Apm(const State& state) const {
    // FFT with zero-padding of the m field
    af::array mfft;
    if (state.Ms_field.isempty()) {
        mfft = af::fftR2C<2>(state.Ms * state.m, af::dim4(nx_expanded(nemesh.nx), ny_expanded(nemesh.ny)));
    } else {
        mfft = af::fftR2C<2>(state.Ms_field * state.m, af::dim4(nx_expanded(nemesh.nx), ny_expanded(nemesh.ny)));
    }

    // Pointwise product
    af::array hfft = af::constant(0.0, nx_expanded(nemesh.nx) / 2 + 1, ny_expanded(nemesh.ny), nemesh.nz, 3, c64);
    for (unsigned i_source = 0; i_source < nemesh.nz; i_source++) {
        for (unsigned i_target = 0; i_target < nemesh.nz; i_target++) {

            int zindex = util::ij2k(i_source, i_target, nemesh.nz);
            af::array nfft;
            if (i_source <= i_target) { // This reflects the data structure of
                                        // newell_nonequi::N_ptr. "<=" choosen
                                        // by definition, could also be ">(=)"
                                        // or just "<" when reconsidering N_ptr
                nfft = Nfft(af::span, af::span, zindex, af::span);
            } else {
                // swap indices to acces symmetric element ij -> ji
                nfft = Nfft(af::span, af::span, zindex, af::span);
                nfft(af::span, af::span, af::span, 2) = -1 * Nfft(af::span, af::span, zindex, 2);
                nfft(af::span, af::span, af::span, 4) = -1 * Nfft(af::span, af::span, zindex, 4);
            }

            hfft(af::span, af::span, i_target, 0) +=
                nfft(af::span, af::span, af::span, 0) * mfft(af::span, af::span, i_source, 0) +
                nfft(af::span, af::span, af::span, 1) * mfft(af::span, af::span, i_source, 1) +
                nfft(af::span, af::span, af::span, 2) * mfft(af::span, af::span, i_source, 2);
            hfft(af::span, af::span, i_target, 1) +=
                nfft(af::span, af::span, af::span, 1) * mfft(af::span, af::span, i_source, 0) +
                nfft(af::span, af::span, af::span, 3) * mfft(af::span, af::span, i_source, 1) +
                nfft(af::span, af::span, af::span, 4) * mfft(af::span, af::span, i_source, 2);
            hfft(af::span, af::span, i_target, 2) +=
                nfft(af::span, af::span, af::span, 2) * mfft(af::span, af::span, i_source, 0) +
                nfft(af::span, af::span, af::span, 4) * mfft(af::span, af::span, i_source, 1) +
                nfft(af::span, af::span, af::span, 5) * mfft(af::span, af::span, i_source, 2);
        }
    }

    af::array one_over_tau_vec = af::array(1, 1, nemesh.nz, 1, f64);
    for (unsigned i = 0; i < nemesh.nz; i++) {
        one_over_tau_vec(0, 0, i, 0) = 1. / (nemesh.dx * nemesh.dy * nemesh.z_spacing[i]);
    }
    one_over_tau_vec = af::tile(one_over_tau_vec, nemesh.nx, nemesh.ny, 1, 3);

    // IFFT reversing padding
    af::array h_field;
    h_field = af::fftC2R<2>(hfft);
    return one_over_tau_vec *
           h_field(af::seq(0, nx_expanded(nemesh.nx) / 2 - 1), af::seq(0, ny_expanded(nemesh.ny) / 2 - 1));
}

namespace newell_nonequi {

double f(double x, double y, double z) {
    x = fabs(x);
    y = fabs(y);
    z = fabs(z);
    const double R = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    const double xx = pow(x, 2);
    const double yy = pow(y, 2);
    const double zz = pow(z, 2);

    double result = 1.0 / 6.0 * (2.0 * xx - yy - zz) * R;
    if (xx + zz > 0)
        result += y / 2.0 * (zz - xx) * asinh(y / (sqrt(xx + zz)));
    if (xx + yy > 0)
        result += z / 2.0 * (yy - xx) * asinh(z / (sqrt(xx + yy)));
    if (x * R > 0)
        result += -x * y * z * atan(y * z / (x * R));
    return result;
}

double g(double x, double y, double z) {
    z = fabs(z);
    const double R = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    const double xx = pow(x, 2);
    const double yy = pow(y, 2);
    const double zz = pow(z, 2);

    double result = -x * y * R / 3.0;
    if (xx + yy > 0)
        result += x * y * z * asinh(z / (sqrt(xx + yy)));
    if (yy + zz > 0)
        result += y / 6.0 * (3.0 * zz - yy) * asinh(x / (sqrt(yy + zz)));
    if (xx + zz > 0)
        result += x / 6.0 * (3.0 * zz - xx) * asinh(y / (sqrt(xx + zz)));
    if (z * R > 0)
        result += -pow(z, 3) / 6.0 * atan(x * y / (z * R));
    if (y * R != 0)
        result += -z * yy / 2.0 * atan(x * z / (y * R));
    if (x * R != 0)
        result += -z * xx / 2.0 * atan(y * z / (x * R));
    return result;
}

// double F2(const double x, const double y, const double z){
//    return f(x, y, z);
//    //Last three terms cancel out//return f(x, y, z) - f(x, 0, z) - f(x, y, 0)
//    + f(x, 0, 0);
//}

double F1(const double x, const double y, const double z, const double dz, const double dZ) {
    return f(x, y, z + dZ) - f(x, y, z) - f(x, y, z - dz + dZ) + f(x, y, z - dz);
}

double F0(const double x, const double y, const double z, const double dy, const double dY, const double dz,
          const double dZ) {
    return F1(x, y + dY, z, dz, dZ) - F1(x, y, z, dz, dZ) - F1(x, y - dy + dY, z, dz, dZ) + F1(x, y - dy, z, dz, dZ);
}

double Nxx(const double x, const double y, const double z, const double dx, const double dy, const double dz,
           const double dX, const double dY, const double dZ) {
    // x, y, z is vector from source cuboid to target cuboid
    // dx, dy, dz are source cuboid dimensions
    // dX, dY, dZ are target cuboid dimensions
    // const double tau = dX * dY * dZ;// Defining dX, dY, dZ as target cuboid
    // (one could alternatively choose dx, dy, dz with implications on x, y, z)
    // return -1./(4.0 * M_PI * tau) * (
    return -1. / (4.0 * M_PI) *
           (F0(x, y, z, dy, dY, dz, dZ) - F0(x - dx, y, z, dy, dY, dz, dZ) - F0(x + dX, y, z, dy, dY, dz, dZ) +
            F0(x - dx + dX, y, z, dy, dY, dz, dZ));
}

// double G2(const double x, const double y, const double z){
//    return g(x, y, z);
//    //return g(x, y, z) - g(x, y, 0);
//    //return g(x, y, z) - g(x, 0, z) - g(x, y, 0) + g(x, 0, 0);
//}

double G1(const double x, const double y, const double z, const double dz, const double dZ) {
    return g(x, y, z + dZ) - g(x, y, z) - g(x, y, z - dz + dZ) + g(x, y, z - dz);
}

double G0(const double x, const double y, const double z, const double dy, const double dY, const double dz,
          const double dZ) {
    return G1(x, y + dY, z, dz, dZ) - G1(x, y, z, dz, dZ) - G1(x, y - dy + dY, z, dz, dZ) + G1(x, y - dy, z, dz, dZ);
}

double Nxy(const double x, const double y, const double z, const double dx, const double dy, const double dz,
           const double dX, const double dY, const double dZ) {
    // x, y, z is vector from source cuboid to target cuboid
    // dx, dy, dz are source cuboid dimensions
    // dX, dY, dZ are target cuboid dimensions
    // const double tau = dX * dY * dZ;// Defining dX, dY, dZ as target cuboid
    // (one could alternatively choose dx, dy, dz with implications on x, y, z)
    // return -1./(4.0 * M_PI * tau) * (
    return -1. / (4.0 * M_PI) *
           (G0(x, y, z, dy, dY, dz, dZ) - G0(x - dx, y, z, dy, dY, dz, dZ) - G0(x + dX, y, z, dy, dY, dz, dZ) +
            G0(x - dx + dX, y, z, dy, dY, dz, dZ));
}

double nonequi_index_distance(const std::vector<double>& spacings, const unsigned i, const unsigned j,
                              const bool verbose) {
    // Calculates the signed distance beween elements by summing up i < j:
    // sum_(k=i)^(j-1)[spacings[k]] or i > j: sum_(k=j)^(i-1)[ - spacings[k]]
    // Note that spacings[spacings.size()-1] is never used
    if (verbose and (i == spacings.size() or j == spacings.size())) {
        printf("%s in nonequi_index_distance: index == vector.size(), the "
               "distance includes the last element which is not wanted "
               "behaviour\n",
               color_string::warning());
    }

    double result = 0;
    if (i > j) {
        for (unsigned k = i; k > j; k--) {
            result -= spacings[k - 1];
        }
    } else {
        for (unsigned k = i; k < j; k++) {
            result += spacings[k];
        }
    }
    return result;
}

void init_N(const NonequiMesh& nemesh, std::vector<double>& N, unsigned ix_start, unsigned ix_end) {
    for (unsigned ix = ix_start; ix < ix_end; ix++) {
        const int jx = (ix + nx_expanded(nemesh.nx) / 2) % nx_expanded(nemesh.nx) - nx_expanded(nemesh.nx) / 2;
        for (unsigned iy = 0; iy < ny_expanded(nemesh.ny); iy++) {
            const int jy = (iy + ny_expanded(nemesh.ny) / 2) % ny_expanded(nemesh.ny) - ny_expanded(nemesh.ny) / 2;
            for (unsigned i_source = 0; i_source < nemesh.nz; i_source++) {
                for (unsigned i_target = 0; i_target < nemesh.nz; i_target++) {

                    if (i_source <= i_target) {
                        const int idx = 6 * (util::ij2k(i_source, i_target, nemesh.nz) +
                                             ((nemesh.nz * (nemesh.nz + 1)) / 2) * (iy + ny_expanded(nemesh.ny) * ix));
                        // std::cout << "idx=" << idx << " of " <<
                        // nx_expanded(nemesh.nx) * ny_expanded(nemesh.ny) * (nemesh.nz *
                        // (nemesh.nz + 1))/2 * 6 << std::endl;
                        const double x = nemesh.dx * (double)jx;
                        const double y = nemesh.dy * (double)jy;
                        const double z = nonequi_index_distance(nemesh.z_spacing, i_source, i_target, true);

                        N[idx + 0] = newell_nonequi::Nxx(x, y, z, nemesh.dx, nemesh.dy, nemesh.z_spacing[i_source],
                                                         nemesh.dx, nemesh.dy, nemesh.z_spacing[i_target]);
                        N[idx + 1] = newell_nonequi::Nxy(x, y, z, nemesh.dx, nemesh.dy, nemesh.z_spacing[i_source],
                                                         nemesh.dx, nemesh.dy, nemesh.z_spacing[i_target]);
                        N[idx + 2] = newell_nonequi::Nxy(x, z, y, nemesh.dx, nemesh.z_spacing[i_source], nemesh.dy,
                                                         nemesh.dx, nemesh.z_spacing[i_target], nemesh.dy);
                        N[idx + 3] = newell_nonequi::Nxx(y, z, x, nemesh.dy, nemesh.z_spacing[i_source], nemesh.dx,
                                                         nemesh.dy, nemesh.z_spacing[i_target], nemesh.dx);
                        N[idx + 4] = newell_nonequi::Nxy(y, z, x, nemesh.dy, nemesh.z_spacing[i_source], nemesh.dx,
                                                         nemesh.dy, nemesh.z_spacing[i_target], nemesh.dx);
                        N[idx + 5] = newell_nonequi::Nxx(z, x, y, nemesh.z_spacing[i_source], nemesh.dx, nemesh.dy,
                                                         nemesh.z_spacing[i_target], nemesh.dx, nemesh.dy);
                    }
                }
            }
        }
    }
}

af::array calculate_N(const NonequiMesh& nemesh, unsigned nthreads) {
    std::vector<double> N_values(nx_expanded(nemesh.nx) * ny_expanded(nemesh.ny) * (nemesh.nz * (nemesh.nz + 1)) / 2 *
                                 6);

    std::vector<std::thread> t;
    for (unsigned i = 0; i < nthreads; i++) {
        unsigned ix_start = i * (double)nx_expanded(nemesh.nx) / nthreads;
        unsigned ix_end = (i + 1) * (double)nx_expanded(nemesh.nx) / nthreads;
        t.push_back(std::thread(init_N, std::ref(nemesh), std::ref(N_values), ix_start, ix_end));
    }

    for (unsigned i = 0; i < nthreads; i++) {
        t[i].join();
    }
    af::array Naf(6, (nemesh.nz * (nemesh.nz + 1)) / 2, ny_expanded(nemesh.ny), nx_expanded(nemesh.nx),
                  N_values.data());
    Naf = af::reorder(Naf, 3, 2, 1, 0);
    Naf = af::fftR2C<2>(Naf);
    return Naf;
}

// wrappered calculate_N
af::array calculate_N(const NonequiMesh& nemesh, unsigned nthreads, bool verbose) {
    af::timer timer = af::timer::start();
    if (verbose) {
        printf("%s Starting Nonequidistant Demag Tensor Assembly on %u out of %u threads.\n", color_string::info(), nthreads,
               std::thread::hardware_concurrency());
    }
    af::array result = calculate_N(nemesh, nthreads);
    if (verbose) {
        printf("time demag init [af-s]: %f\n", af::timer::stop(timer));
    }
    return result;
}

} // namespace newell_nonequi

namespace nonequi_util {
std::string to_string(const NonequiMesh& nemesh) {
    std::string dz_string;
    for (auto const& dz : nemesh.z_spacing) {
        dz_string.append(std::to_string(1e9 * dz));
    }
    return "n0exp_" + std::to_string(nx_expanded(nemesh.nx)) + "_n1exp_" + std::to_string(ny_expanded(nemesh.ny)) +
           "_nz_" + std::to_string(nemesh.nz) + "_dx_" + std::to_string(1e9 * nemesh.dx) + "_dy_" +
           std::to_string(1e9 * nemesh.dy) + "_dz_" + dz_string;
}

void warn_if_maxprime_lt_13(unsigned n, const std::string& ni) {
    if (util::max_of_prime_factors(n) > 13) {
        std::cout << color_string::warning() << " NonequiDemagField::NonequiDemagField: maximum prime factor of nemesh." << ni << "="
                  << n << " is " << util::max_of_prime_factors(n)
                  << ", which is > 13. FFT on the OpenCL backend only supports dimensions with the maximum prime "
                     "factor <= 13. Please use either the CUDA or CPU backend or choose an alternative discretization "
                     "where max_prime(n) <= 13."
                  << std::endl;
    }
}
void warn_if_maxprime_lt_13(const NonequiMesh& nemesh) {
    if (af::getActiveBackend() == AF_BACKEND_OPENCL) {
        warn_if_maxprime_lt_13(nemesh.nx, "nx");
        warn_if_maxprime_lt_13(nemesh.ny, "ny");
        warn_if_maxprime_lt_13(nemesh.nz, "nz");
    }
}
} // namespace nonequi_util

af::array get_Nfft(const NonequiMesh& nemesh, bool verbose, bool caching, unsigned nthreads) {
    nonequi_util::warn_if_maxprime_lt_13(nemesh);

    if (caching == false) {
        return newell_nonequi::calculate_N(nemesh, nthreads, verbose);
    } else {
        util::CacheManager cm{verbose};
        const std::string nfft_id = nonequi_util::to_string(nemesh);
        auto optional_Nfft = cm.get_array_if_existent(nfft_id);
        if (optional_Nfft) {
            return optional_Nfft.value();
        } else {
            auto result = newell_nonequi::calculate_N(nemesh, nthreads, verbose);
            cm.write_array(result, nfft_id);
            return result;
        }
    }
}

NonequiDemagField::NonequiDemagField(const NonequiMesh& nemesh, bool verbose, bool caching, unsigned nthreads)
    : NonequiTerm(nemesh),
      Nfft(get_Nfft(nemesh, verbose, caching, nthreads > 0 ? nthreads : std::thread::hardware_concurrency())) {}

} // namespace magnumafcpp
