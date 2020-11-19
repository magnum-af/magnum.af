#include "field_terms/micro_demag.hpp"
#include "func.hpp"
#include "misc.hpp"
#include "util/af_overloads.hpp"
#include "util/cache_manager.hpp"
#include "util/prime_factors.hpp"
#include <string>
#include <thread>

namespace magnumafcpp {

void DemagField::print_Nfft() const { af::print("Nfft=", Nfft); }

af::array calculate_N(int n0_exp, int n1_exp, int n2_exp, double dx, double dy, double dz, unsigned nthreads);
af::array calculate_N(Mesh mesh, bool verbose, unsigned nthreads) {
    af::timer demagtimer = af::timer::start();
    if (verbose) {
        printf("%s Starting Demag Tensor Assembly on %u out of %u threads.\n", Info(), nthreads,
               std::thread::hardware_concurrency());
    }
    auto result = calculate_N(mesh.n0_exp, mesh.n1_exp, mesh.n2_exp, mesh.dx, mesh.dy, mesh.dz, nthreads);
    if (verbose) {
        printf("%s Initialized demag tensor in %f [af-s]\n", Info(), af::timer::stop(demagtimer));
    }
    return result;
}

void warn_if_maxprime_lt_13(unsigned n, std::string ni) {
    if (util::max_of_prime_factors(n) > 13) {
        std::cout << Warning() << " DemagField::DemagField: maximum prime factor of mesh." << ni << "=" << n << " is "
                  << util::max_of_prime_factors(n)
                  << ", which is > 13. FFT on the OpenCL backend only supports dimensions with the maximum prime "
                     "factor <= 13. Please choose an alternative discretization where max_prime(n) <= 13."
                  << std::endl;
    }
}
void warn_if_maxprime_lt_13(Mesh mesh) {
    warn_if_maxprime_lt_13(mesh.n0, "nx");
    warn_if_maxprime_lt_13(mesh.n1, "ny");
    warn_if_maxprime_lt_13(mesh.n2, "nz");
}

std::string to_string(const Mesh& mesh) {
    return "Nfft_n0exp_" + std::to_string(mesh.n0_exp) + "_n1exp_" + std::to_string(mesh.n1_exp) + "_n2exp_" +
           std::to_string(mesh.n2_exp) + "_dx_" + std::to_string(1e9 * mesh.dx) + "_dy_" +
           std::to_string(1e9 * mesh.dy) + "_dz_" + std::to_string(1e9 * mesh.dz);
}

af::array get_Nfft(Mesh mesh, bool verbose, bool caching, unsigned nthreads) {
    if (af::getActiveBackend() == AF_BACKEND_OPENCL) {
        warn_if_maxprime_lt_13(mesh);
    }

    if (caching == false) {
        return calculate_N(mesh, verbose, nthreads);
    } else {
        util::CacheManager cm{verbose};
        const std::string nfft_id = to_string(mesh);
        auto optional_Nfft = cm.get_array_if_existent(nfft_id);
        if (optional_Nfft) {
            return optional_Nfft.value();
        } else {
            auto result = calculate_N(mesh, verbose, nthreads);
            cm.write_array(result, nfft_id);
            return result;
        }
    }
}

DemagField::DemagField(Mesh mesh, bool verbose, bool caching, unsigned in_nthreads)
    : Nfft(::magnumafcpp::get_Nfft(mesh, verbose, caching,
                                   in_nthreads > 0 ? in_nthreads : std::thread::hardware_concurrency())) {}

af::array DemagField::h(const State& state) {
    af::timer timer_demagsolve = af::timer::start();

    // Converting Nfft from c64 to c32 once if state.m.type() == f32
    if (Nfft.type() == af::dtype::c64 and state.m.type() == af::dtype::f32) {
        std::cout << "DemagField::h: state.m is of type " << state.m.type() << ", converting Nfft type from "
                  << Nfft.type() << " to " << af::dtype::c32 << std::endl;
        Nfft = Nfft.as(af::dtype::c32);
    }

    // FFT with zero-padding of the m field
    af::array mfft;
    if (state.mesh.n2_exp == 1) {
        if (state.Ms_field.isempty())
            mfft = af::fftR2C<2>(state.Ms * state.m, af::dim4(state.mesh.n0_exp, state.mesh.n1_exp));
        else
            mfft = af::fftR2C<2>(state.Ms_field * state.m, af::dim4(state.mesh.n0_exp, state.mesh.n1_exp));
    } else {
        if (state.Ms_field.isempty())
            mfft = af::fftR2C<3>(state.Ms * state.m, af::dim4(state.mesh.n0_exp, state.mesh.n1_exp, state.mesh.n2_exp));
        else
            mfft = af::fftR2C<3>(state.Ms_field * state.m,
                                 af::dim4(state.mesh.n0_exp, state.mesh.n1_exp, state.mesh.n2_exp));
    }

    // Pointwise product
    af::array hfft = af::array(state.mesh.n0_exp / 2 + 1, state.mesh.n1_exp, state.mesh.n2_exp, 3, Nfft.type());
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

    // IFFT reversing padding
    af::array h_field;
    if (state.mesh.n2_exp == 1) {
        h_field = af::fftC2R<2>(hfft);
        if (state.afsync)
            af::sync();
        cpu_time += af::timer::stop(timer_demagsolve);
        return h_field(af::seq(0, state.mesh.n0_exp / 2 - 1), af::seq(0, state.mesh.n1_exp / 2 - 1));
    } else {
        h_field = af::fftC2R<3>(hfft);
        if (state.afsync)
            af::sync();
        cpu_time += af::timer::stop(timer_demagsolve);
        return h_field(af::seq(0, state.mesh.n0_exp / 2 - 1), af::seq(0, state.mesh.n1_exp / 2 - 1),
                       af::seq(0, state.mesh.n2_exp / 2 - 1), af::span);
    }
}

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
    if (xx + zz > 0)
        result += y / 2.0 * (zz - xx) * asinh(y / (sqrt(xx + zz)));
    if (xx + yy > 0)
        result += z / 2.0 * (yy - xx) * asinh(z / (sqrt(xx + yy)));
    if (x * R > 0)
        result += -x * y * z * atan(y * z / (x * R));
    return result;
}

double g(const double x, const double y, double z) {
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

struct LoopInfo {
    LoopInfo() {}
    LoopInfo(int i0_start, int i0_end, int n0_exp, int n1_exp, int n2_exp, double dx, double dy, double dz)
        : i0_start(i0_start), i0_end(i0_end), n0_exp(n0_exp), n1_exp(n1_exp), n2_exp(n2_exp), dx(dx), dy(dy), dz(dz) {}
    int i0_start;
    int i0_end;
    int n0_exp;
    int n1_exp;
    int n2_exp;
    double dx;
    double dy;
    double dz;
};

void setup_N(const LoopInfo& loopinfo, std::vector<double>& N) {
    for (int i0 = loopinfo.i0_start; i0 < loopinfo.i0_end; i0++) {
        const int j0 = (i0 + loopinfo.n0_exp / 2) % loopinfo.n0_exp - loopinfo.n0_exp / 2;
        for (int i1 = 0; i1 < loopinfo.n1_exp; i1++) {
            const int j1 = (i1 + loopinfo.n1_exp / 2) % loopinfo.n1_exp - loopinfo.n1_exp / 2;
            for (int i2 = 0; i2 < loopinfo.n2_exp; i2++) {
                const int j2 = (i2 + loopinfo.n2_exp / 2) % loopinfo.n2_exp - loopinfo.n2_exp / 2;
                const int idx = 6 * (i2 + loopinfo.n2_exp * (i1 + loopinfo.n1_exp * i0));
                N[idx + 0] = newell::Nxx(j0, j1, j2, loopinfo.dx, loopinfo.dy, loopinfo.dz);
                N[idx + 1] = newell::Nxy(j0, j1, j2, loopinfo.dx, loopinfo.dy, loopinfo.dz);
                N[idx + 2] = newell::Nxy(j0, j2, j1, loopinfo.dx, loopinfo.dz, loopinfo.dy);
                N[idx + 3] = newell::Nxx(j1, j2, j0, loopinfo.dy, loopinfo.dz, loopinfo.dx);
                N[idx + 4] = newell::Nxy(j1, j2, j0, loopinfo.dy, loopinfo.dz, loopinfo.dx);
                N[idx + 5] = newell::Nxx(j2, j0, j1, loopinfo.dz, loopinfo.dx, loopinfo.dy);
            }
        }
    }
}
} // namespace newell

af::array calculate_N(int n0_exp, int n1_exp, int n2_exp, double dx, double dy, double dz, unsigned nthreads) {
    std::vector<newell::LoopInfo> loopinfo;
    for (unsigned i = 0; i < nthreads; i++) {
        unsigned start = i * (double)n0_exp / nthreads;
        unsigned end = (i + 1) * (double)n0_exp / nthreads;
        loopinfo.push_back(newell::LoopInfo(start, end, n0_exp, n1_exp, n2_exp, dx, dy, dz));
    }

    std::vector<std::thread> t;
    std::vector<double> N_values(n0_exp * n1_exp * n2_exp * 6);
    for (unsigned i = 0; i < nthreads; i++) {
        t.push_back(std::thread(newell::setup_N, std::ref(loopinfo[i]), std::ref(N_values)));
    }

    for (unsigned i = 0; i < nthreads; i++) {
        t[i].join();
    }

    af::array Naf(6, n2_exp, n1_exp, n0_exp, N_values.data());
    Naf = af::reorder(Naf, 3, 2, 1, 0);

    if (n2_exp == 1) {
        Naf = af::fftR2C<2>(Naf);
    } else {
        Naf = af::fftR2C<3>(Naf);
    }
    return Naf;
}
} // namespace magnumafcpp
