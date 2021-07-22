#include "field_terms/micro/demag_field_pbc.hpp"
#include "util/af_overloads.hpp"
#include "util/color_string.hpp"
#include "util/prime_factors.hpp"
#include <cmath>

namespace magnumaf {

namespace {
af::array calc_m_mfft(const State& state) {
    if (state.mesh.nz == 1) {
        if (state.Ms_field.isempty()) {
            return af::fftR2C<2>(state.Ms * state.m);
        } else {
            return af::fftR2C<2>(state.Ms_field * state.m);
        }
    } else {
        if (state.Ms_field.isempty()) {
            return af::fftR2C<3>(state.Ms * state.m);
        } else {
            return af::fftR2C<3>(state.Ms_field * state.m);

}
    }
}

auto fftC2R_dim2switch = [](const af::array& h_fft) {
    if (h_fft.dims(2) == 1) {
        return af::fftC2R<2>(h_fft);
    } else {
        return af::fftC2R<3>(h_fft);
    }
};

void warn_if_maxprime_gt_13_opencl(std::size_t n, const std::string& ni) {
    if (util::max_of_prime_factors(n) > 13) {
        std::cout << color_string::warning() << "DemagFieldPBC: maximum prime factor of 'mesh." << ni << "'=" << n
                  << " is " << util::max_of_prime_factors(n)
                  << ", which is > 13. FFT on the OpenCL backend only supports dimensions with the maximum prime "
                     "factor <= 13. Please use either the CUDA or CPU backend or choose an alternative discretization "
                     "where max_prime(n) <= 13."
                  << std::endl;
    }
}

void warn_if_maxprime_gt_13_opencl(const Mesh& mesh) {
    if (af::getActiveBackend() == AF_BACKEND_OPENCL) {
        warn_if_maxprime_gt_13_opencl(mesh.nx, "nx");
        warn_if_maxprime_gt_13_opencl(mesh.ny, "ny");
        warn_if_maxprime_gt_13_opencl(mesh.nz, "nz");
    }
}

} // namespace

af::array DemagFieldPBC::impl_H_in_Apm(const State& state) const {
    warn_if_maxprime_gt_13_opencl(state.mesh);
    const auto m_fft = calc_m_mfft(state);
    const auto nx_full = state.mesh.nx;
    const auto nx = m_fft.dims(0); // nx = nx_full / 2 + 1
    const auto ny = m_fft.dims(1);
    const auto nz = m_fft.dims(2);
    const auto dx = state.mesh.dx;
    const auto dy = state.mesh.dy;
    const auto dz = state.mesh.dz;

    const auto kx = (2. * af::Pi * af::range(nx) / nx_full);
    const auto ky = af::reorder((2. * af::Pi * af::range(ny) / ny), 1, 0);
    const auto kz = af::reorder((2. * af::Pi * af::range(nz) / nz), 1, 2, 0);
    const auto kx_ = [kx, ny, nz] { return af::tile(kx, 1, ny, nz); };
    const auto ky_ = [ky, nx, nz] { return af::tile(ky, nx, 1, nz); };
    const auto kz_ = [kz, nx, ny] { return af::tile(kz, nx, ny, 1); };

    // Make pure complex array, NOTE: af::complex(a) returns complex (a, 0), not (0, a)
    const auto to_j = [](const af::array& a) { return af::complex(af::constant(0., a.dims(), a.type()), a); };

    const auto div_fft = (1. - af::exp(-1. * to_j(kx_()))) * m_fft(af::span, af::span, af::span, 0) / dx +
                         (1. - af::exp(-1. * to_j(ky_()))) * m_fft(af::span, af::span, af::span, 1) / dy +
                         (1. - af::exp(-1. * to_j(kz_()))) * m_fft(af::span, af::span, af::span, 2) / dz;

    auto u_fft = -div_fft / (4. / std::pow(dx, 2) * af::pow(af::sin(kx_() / 2.), 2) +
                             4. / std::pow(dy, 2) * af::pow(af::sin(ky_() / 2.), 2) +
                             4. / std::pow(dz, 2) * af::pow(af::sin(kz_() / 2.), 2));
    u_fft(0, 0, 0) = 0; // self-interaction term has division by zero, we overwrite nans with zero

    const auto h_fft_x = (1. - af::exp(to_j(kx_()))) * u_fft / dx;
    const auto h_fft_y = (1. - af::exp(to_j(ky_()))) * u_fft / dy;
    const auto h_fft_z = (1. - af::exp(to_j(kz_()))) * u_fft / dz;
    const auto h_fft = af::join(3, h_fft_x, h_fft_y, h_fft_z);
    return fftC2R_dim2switch(h_fft);
}
} // namespace magnumaf
