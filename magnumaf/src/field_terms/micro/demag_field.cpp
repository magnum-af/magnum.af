#include "field_terms/micro/demag_field.hpp"
#include "util/af_overloads.hpp"
#include "util/cache_manager.hpp"
#include "util/color_string.hpp"
#include "util/prime_factors.hpp"
#include "util/util.hpp"
#include <string>
#include <thread>

namespace magnumaf {

namespace {
// Expanded cell sizes for demag FFT
unsigned nx_exp(unsigned nx) { return 2 * nx; }
unsigned ny_exp(unsigned ny) { return 2 * ny; }
unsigned nz_exp(unsigned nz) { return (nz == 1) ? 1 : 2 * nz; }
} // namespace

namespace newell {
double f(double x, double y, double z) {
    x = fabs(x);
    y = fabs(y);
    z = fabs(z);
    const double R = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    const double xx = pow(x, 2);
    const double yy = pow(y, 2);
    const double zz = pow(z, 2);

    double result = 1.0 / 6.0 * (2.0 * xx - yy - zz) * R;
    if (xx + zz > 0) {
        result += y / 2.0 * (zz - xx) * asinh(y / (sqrt(xx + zz)));
    }
    if (xx + yy > 0) {
        result += z / 2.0 * (yy - xx) * asinh(z / (sqrt(xx + yy)));
    }
    if (x * R > 0) {
        result += -x * y * z * atan(y * z / (x * R));
    }
    return result;
}

double g(const double x, const double y, double z) {
    z = fabs(z);
    const double R = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    const double xx = pow(x, 2);
    const double yy = pow(y, 2);
    const double zz = pow(z, 2);

    double result = -x * y * R / 3.0;
    if (xx + yy > 0) {
        result += x * y * z * asinh(z / (sqrt(xx + yy)));
    }
    if (yy + zz > 0) {
        result += y / 6.0 * (3.0 * zz - yy) * asinh(x / (sqrt(yy + zz)));
    }
    if (xx + zz > 0) {
        result += x / 6.0 * (3.0 * zz - xx) * asinh(y / (sqrt(xx + zz)));
    }
    if (z * R > 0) {
        result += -pow(z, 3) / 6.0 * atan(x * y / (z * R));
    }
    if (y * R != 0) {
        result += -z * yy / 2.0 * atan(x * z / (y * R));
    }
    if (x * R != 0) {
        result += -z * xx / 2.0 * atan(y * z / (x * R));
    }
    return result;
}

double Nxx(const int ix, const int iy, const int iz, const double dx, const double dy, const double dz) {
    const double x = dx * ix;
    const double y = dy * iy;
    const double z = dz * iz;
    double result = 8.0 * f(x, y, z) - 4.0 * f(x + dx, y, z) - 4.0 * f(x - dx, y, z) - 4.0 * f(x, y + dy, z) -
                    4.0 * f(x, y - dy, z) - 4.0 * f(x, y, z + dz) - 4.0 * f(x, y, z - dz) + 2.0 * f(x + dx, y + dy, z) +
                    2.0 * f(x + dx, y - dy, z) + 2.0 * f(x - dx, y + dy, z) + 2.0 * f(x - dx, y - dy, z) +
                    2.0 * f(x + dx, y, z + dz) + 2.0 * f(x + dx, y, z - dz) + 2.0 * f(x - dx, y, z + dz) +
                    2.0 * f(x - dx, y, z - dz) + 2.0 * f(x, y + dy, z + dz) + 2.0 * f(x, y + dy, z - dz) +
                    2.0 * f(x, y - dy, z + dz) + 2.0 * f(x, y - dy, z - dz) - 1.0 * f(x + dx, y + dy, z + dz) -
                    1.0 * f(x + dx, y + dy, z - dz) - 1.0 * f(x + dx, y - dy, z + dz) -
                    1.0 * f(x + dx, y - dy, z - dz) - 1.0 * f(x - dx, y + dy, z + dz) -
                    1.0 * f(x - dx, y + dy, z - dz) - 1.0 * f(x - dx, y - dy, z + dz) - 1.0 * f(x - dx, y - dy, z - dz);
    return -result / (4.0 * M_PI * dx * dy * dz);
}

double Nxy(const int ix, const int iy, const int iz, const double dx, const double dy, const double dz) {
    const double x = dx * ix;
    const double y = dy * iy;
    const double z = dz * iz;
    double result = 8.0 * g(x, y, z) - 4.0 * g(x + dx, y, z) - 4.0 * g(x - dx, y, z) - 4.0 * g(x, y + dy, z) -
                    4.0 * g(x, y - dy, z) - 4.0 * g(x, y, z + dz) - 4.0 * g(x, y, z - dz) + 2.0 * g(x + dx, y + dy, z) +
                    2.0 * g(x + dx, y - dy, z) + 2.0 * g(x - dx, y + dy, z) + 2.0 * g(x - dx, y - dy, z) +
                    2.0 * g(x + dx, y, z + dz) + 2.0 * g(x + dx, y, z - dz) + 2.0 * g(x - dx, y, z + dz) +
                    2.0 * g(x - dx, y, z - dz) + 2.0 * g(x, y + dy, z + dz) + 2.0 * g(x, y + dy, z - dz) +
                    2.0 * g(x, y - dy, z + dz) + 2.0 * g(x, y - dy, z - dz) - 1.0 * g(x + dx, y + dy, z + dz) -
                    1.0 * g(x + dx, y + dy, z - dz) - 1.0 * g(x + dx, y - dy, z + dz) -
                    1.0 * g(x + dx, y - dy, z - dz) - 1.0 * g(x - dx, y + dy, z + dz) -
                    1.0 * g(x - dx, y + dy, z - dz) - 1.0 * g(x - dx, y - dy, z + dz) - 1.0 * g(x - dx, y - dy, z - dz);
    result = -result / (4.0 * M_PI * dx * dy * dz);
    return result;
}

void setup_N(const Mesh& mesh, std::vector<double>& N, unsigned ix_start, unsigned ix_end) {
    for (unsigned i0 = ix_start; i0 < ix_end; i0++) {
        const int j0 = (i0 + nx_exp(mesh.nx) / 2) % nx_exp(mesh.nx) - nx_exp(mesh.nx) / 2;
        for (unsigned i1 = 0; i1 < ny_exp(mesh.ny); i1++) {
            const int j1 = (i1 + ny_exp(mesh.ny) / 2) % ny_exp(mesh.ny) - ny_exp(mesh.ny) / 2;
            for (unsigned i2 = 0; i2 < nz_exp(mesh.nz); i2++) {
                const int j2 = (i2 + nz_exp(mesh.nz) / 2) % nz_exp(mesh.nz) - nz_exp(mesh.nz) / 2;
                const int idx = 6 * (i2 + nz_exp(mesh.nz) * (i1 + ny_exp(mesh.ny) * i0));
                N[idx + 0] = newell::Nxx(j0, j1, j2, mesh.dx, mesh.dy, mesh.dz);
                N[idx + 1] = newell::Nxy(j0, j1, j2, mesh.dx, mesh.dy, mesh.dz);
                N[idx + 2] = newell::Nxy(j0, j2, j1, mesh.dx, mesh.dz, mesh.dy);
                N[idx + 3] = newell::Nxx(j1, j2, j0, mesh.dy, mesh.dz, mesh.dx);
                N[idx + 4] = newell::Nxy(j1, j2, j0, mesh.dy, mesh.dz, mesh.dx);
                N[idx + 5] = newell::Nxx(j2, j0, j1, mesh.dz, mesh.dx, mesh.dy);
            }
        }
    }
}
} // namespace newell

namespace {
af::array calculate_N_real(const Mesh& mesh, unsigned nthreads) {
    std::vector<double> N_values(nx_exp(mesh.nx) * ny_exp(mesh.ny) * nz_exp(mesh.nz) * 6);
    std::vector<std::thread> t;

    for (unsigned i = 0; i < nthreads; i++) {
        unsigned ix_start = i * static_cast<double>(nx_exp(mesh.nx)) / nthreads;
        unsigned ix_end = (i + 1) * static_cast<double>(nx_exp(mesh.nx)) / nthreads;
        t.emplace_back(newell::setup_N, std::ref(mesh), std::ref(N_values), ix_start, ix_end);
    }

    for (unsigned i = 0; i < nthreads; i++) {
        t[i].join();
    }

    af::array Naf(6, nz_exp(mesh.nz), ny_exp(mesh.ny), nx_exp(mesh.nx), N_values.data());
    Naf = af::reorder(Naf, 3, 2, 1, 0);
    return Naf;
}

af::array calculate_Nfft(const Mesh& mesh, unsigned nthreads) {
    const af::array Naf = calculate_N_real(mesh, nthreads);
    return nz_exp(mesh.nz) == 1 ? af::fftR2C<2>(Naf) : af::fftR2C<3>(Naf);
}

af::array calculate_Nfft(Mesh mesh, bool verbose, unsigned nthreads) {
    af::timer demagtimer = af::timer::start();
    if (verbose) {
        printf("%s Starting Demag Tensor Assembly on %u out of %u threads.\n", color_string::info(), nthreads,
               std::thread::hardware_concurrency());
    }
    auto result = calculate_Nfft(mesh, nthreads);
    if (verbose) {
        printf("%s Initialized demag tensor in %f [af-s]\n", color_string::info(), af::timer::stop(demagtimer));
    }
    return result;
}

void warn_if_maxprime_gt_13(std::size_t n, const std::string& ni) {
    if (util::max_of_prime_factors(n) > 13) {
        std::cout << color_string::warning() << "DemagField::DemagField: maximum prime factor of mesh." << ni << "="
                  << n << " is " << util::max_of_prime_factors(n)
                  << ", which is > 13. FFT on the OpenCL backend only supports dimensions with the maximum prime "
                     "factor <= 13. Please use either the CUDA or CPU backend or choose an alternative discretization "
                     "where max_prime(n) <= 13."
                  << std::endl;
    }
}

void warn_if_maxprime_gt_13(const Mesh& mesh) {
    warn_if_maxprime_gt_13(mesh.nx, "nx");
    warn_if_maxprime_gt_13(mesh.ny, "ny");
    warn_if_maxprime_gt_13(mesh.nz, "nz");
}

std::string to_string(const Mesh& mesh) {
    return "Nfft_n0exp_" + std::to_string(nx_exp(mesh.nx)) + "_n1exp_" + std::to_string(ny_exp(mesh.ny)) + "_n2exp_" +
           std::to_string(nz_exp(mesh.nz)) + "_dx_" + std::to_string(1e9 * mesh.dx) + "_dy_" +
           std::to_string(1e9 * mesh.dy) + "_dz_" + std::to_string(1e9 * mesh.dz);
}

af::array get_Nfft(Mesh mesh, bool verbose, bool caching, std::size_t nthreads) {
    nthreads = nthreads > 0 ? nthreads : std::thread::hardware_concurrency();
    if (af::getActiveBackend() == AF_BACKEND_OPENCL) {
        warn_if_maxprime_gt_13(mesh);
    }

    if (caching == false) {
        return calculate_Nfft(mesh, verbose, nthreads);
    } else {
        util::CacheManager cm{verbose};
        const std::string nfft_id = to_string(mesh);
        auto optional_Nfft = cm.get_array_if_existent(nfft_id);
        if (optional_Nfft) {
            return optional_Nfft.value();
        } else {
            auto result = calculate_Nfft(mesh, verbose, nthreads);
            cm.write_array(result, nfft_id);
            return result;
        }
    }
}
} // namespace

af::array DemagField::calc_N(Mesh mesh) const {
    return ::magnumaf::calculate_N_real(mesh, std::thread::hardware_concurrency());
}

DemagField::DemagField(Mesh mesh, bool verbose, bool caching, unsigned nthreads)
    : Nfft(::magnumaf::get_Nfft(mesh, verbose, caching, nthreads)) {}

af::array DemagField::impl_H_in_Apm(const State& state) const {
    // Converting Nfft from c64 to c32 once if state.m.type() == f32
    if (Nfft.type() == af::dtype::c64 and state.m.type() == af::dtype::f32) {
        std::cout << "DemagField::h: state.m is of type " << state.m.type() << ", converting Nfft type from "
                  << Nfft.type() << " to " << af::dtype::c32 << std::endl;
        Nfft = Nfft.as(af::dtype::c32);
    }

    const af::array Ms_times_m =
        state.Ms_field.isempty() ? state.Ms * state.m : state.get_Ms_field_in_vec_dims() * state.m;

    // FFT with zero-padding of the m field
    const af::array mfft =
        nz_exp(state.mesh.nz) == 1
            ? af::fftR2C<2>(Ms_times_m, af::dim4(nx_exp(state.mesh.nx), ny_exp(state.mesh.ny)))
            : af::fftR2C<3>(Ms_times_m, af::dim4(nx_exp(state.mesh.nx), ny_exp(state.mesh.ny), nz_exp(state.mesh.nz)));

    // Pointwise product
    af::array hfft =
        af::array(nx_exp(state.mesh.nx) / 2 + 1, ny_exp(state.mesh.ny), nz_exp(state.mesh.nz), 3, Nfft.type());
    hfft(af::span, af::span, af::span, 0) =
        Nfft(af::span, af::span, af::span, 0) * mfft(af::span, af::span, af::span, 0) +
        Nfft(af::span, af::span, af::span, 1) * mfft(af::span, af::span, af::span, 1) +
        Nfft(af::span, af::span, af::span, 2) * mfft(af::span, af::span, af::span, 2);
    hfft(af::span, af::span, af::span, 1) =
        Nfft(af::span, af::span, af::span, 1) * mfft(af::span, af::span, af::span, 0) +
        Nfft(af::span, af::span, af::span, 3) * mfft(af::span, af::span, af::span, 1) +
        Nfft(af::span, af::span, af::span, 4) * mfft(af::span, af::span, af::span, 2);
    hfft(af::span, af::span, af::span, 2) =
        Nfft(af::span, af::span, af::span, 2) * mfft(af::span, af::span, af::span, 0) +
        Nfft(af::span, af::span, af::span, 4) * mfft(af::span, af::span, af::span, 1) +
        Nfft(af::span, af::span, af::span, 5) * mfft(af::span, af::span, af::span, 2);

    constexpr double legacy_fft_norm_factor = 0.0; // legacy defualt value pre v3.8.1
    // NOTE: using legacy default norm factor of 0.0 instead of current default 1.0 prevents a bug introduced via
    // arrayfire v3.8.1 https://github.com/arrayfire/arrayfire/pull/3178

    // IFFT reversing padding
    const af::array h_field_padded = nz_exp(state.mesh.nz) == 1 ? af::fftC2R<2>(hfft, false, legacy_fft_norm_factor)
                                                                : af::fftC2R<3>(hfft, false, legacy_fft_norm_factor);
    const af::array h_field =
        nz_exp(state.mesh.nz) == 1
            ? h_field_padded(af::seq(0, nx_exp(state.mesh.nx) / 2 - 1), af::seq(0, ny_exp(state.mesh.ny) / 2 - 1))
            : h_field_padded(af::seq(0, nx_exp(state.mesh.nx) / 2 - 1), af::seq(0, ny_exp(state.mesh.ny) / 2 - 1),
                             af::seq(0, nz_exp(state.mesh.nz) / 2 - 1), af::span);
    return h_field;
}

} // namespace magnumaf
