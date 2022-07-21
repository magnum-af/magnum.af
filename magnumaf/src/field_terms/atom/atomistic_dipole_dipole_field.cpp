#include "field_terms/atom/atomistic_dipole_dipole_field.hpp"
#include "util/util.hpp"
#include <vector>

namespace magnumaf {

// Expanded cell sizes for demag FFT
inline unsigned nx_exp(unsigned nx) { return 2 * nx; }
inline unsigned ny_exp(unsigned ny) { return 2 * ny; }
inline unsigned nz_exp(unsigned nz) { return (nz == 1) ? 1 : 2 * nz; }

af::array N_atomistic(int n0_exp, int n1_exp, int n2_exp, double dx, double dy, double dz);

AtomisticDipoleDipoleField::AtomisticDipoleDipoleField(Mesh mesh)
    : Nfft(N_atomistic(nx_exp(mesh.nx), ny_exp(mesh.ny), nz_exp(mesh.nz), mesh.dx, mesh.dy, mesh.dz)) {}

af::array AtomisticDipoleDipoleField::impl_H_in_Apm(const State& state) const {
    // FFT with zero-padding of the m field
    af::array mfft;
    if (nz_exp(state.mesh.nz) == 1) {
        mfft = af::fftR2C<2>(state.m, af::dim4(nx_exp(state.mesh.nx), ny_exp(state.mesh.ny)));
    } else {
        mfft = af::fftR2C<3>(state.m, af::dim4(nx_exp(state.mesh.nx), ny_exp(state.mesh.ny), nz_exp(state.mesh.nz)));
    }

    af::array hfft = af::array(nx_exp(state.mesh.nx) / 2 + 1, ny_exp(state.mesh.ny), nz_exp(state.mesh.nz), 3, c64);
    // Pointwise product
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
    if (nz_exp(state.mesh.nz) == 1) {
        h_field = af::fftC2R<2>(hfft, false, 0.0);
        return state.Ms * h_field(af::seq(0, nx_exp(state.mesh.nx) / 2 - 1),
                                  af::seq(0, ny_exp(state.mesh.ny) / 2 - 1)); // TODO consider p density, then we
                                                                              // have to multip at m before fft
    } else {
        h_field = af::fftC2R<3>(hfft, false, 0.0);
        return state.Ms * h_field(af::seq(0, nx_exp(state.mesh.nx) / 2 - 1), af::seq(0, ny_exp(state.mesh.ny) / 2 - 1),
                                  af::seq(0, nz_exp(state.mesh.nz) / 2 - 1), af::span);
    }
}

af::array N_atomistic(int n0_exp, int n1_exp, int n2_exp, double dx, double dy, double dz) {
    std::vector<double> N(n0_exp * n1_exp * n2_exp * 6);
    // Experimental
    for (int i0 = 0; i0 < n0_exp; i0++) {
        const int j0 = (i0 + n0_exp / 2) % n0_exp - n0_exp / 2;
        for (int i1 = 0; i1 < n1_exp; i1++) {
            const int j1 = (i1 + n1_exp / 2) % n1_exp - n1_exp / 2;
            for (int i2 = 0; i2 < n2_exp; i2++) {
                const int j2 = (i2 + n2_exp / 2) % n2_exp - n2_exp / 2;
                const int idx = 6 * (i2 + n2_exp * (i1 + n1_exp * i0));
                const double rx = j0 * dx;
                const double ry = j1 * dy;
                const double rz = j2 * dz;
                const double r = sqrt(pow(rx, 2) + pow(ry, 2) + pow(rz, 2));
                if (r == 0.) { // TODO repsace with if (j0 == 0 && j1 == 0 && j2
                               // == 0)
                    // std::cout<<"In AtomisticDipoleDipoleField::N_atomistic:
                    // r==0"<<std::endl; std::cout<<"In
                    // AtomisticDipoleDipoleField::setting n to
                    // 1/3."<<std::endl; Accounting for self-interaction (would
                    // be inf, when approximated with sphere -1/3 in diag
                    // TODO check
                    N[idx + 0] = 0.;
                    N[idx + 1] = 0.;
                    N[idx + 2] = 0.;
                    N[idx + 3] = 0.;
                    N[idx + 4] = 0.;
                    N[idx + 5] = 0.;

                    // N[idx+0] = -1./3.;
                    // N[idx+1] = 0.;
                    // N[idx+2] = 0.;
                    // N[idx+3] = -1./3.;
                    // N[idx+4] = 0.;
                    // N[idx+5] = -1./3.;
                } else {
                    // N[idx+0] = 1./(4.*M_PI)*(3.*rx*rx/pow(r, 5) - 1./pow(r,
                    // 3)); N[idx+1] = 1./(4.*M_PI)*(3.*rx*ry/pow(r, 5) );
                    // N[idx+2] = 1./(4.*M_PI)*(3.*rx*rz/pow(r, 5) ); N[idx+3]
                    // = 1./(4.*M_PI)*(3.*ry*ry/pow(r, 5) - 1./pow(r, 3));
                    // N[idx+4] = 1./(4.*M_PI)*(3.*ry*rz/pow(r, 5) ); N[idx+5]
                    // = 1./(4.*M_PI)*(3.*rz*rz/pow(r, 5) - 1./pow(r, 3));
                    N[idx + 0] = 3. * rx * rx / pow(r, 5) - 1. / pow(r, 3);
                    N[idx + 1] = 3. * rx * ry / pow(r, 5);
                    N[idx + 2] = 3. * rx * rz / pow(r, 5);
                    N[idx + 3] = 3. * ry * ry / pow(r, 5) - 1. / pow(r, 3);
                    N[idx + 4] = 3. * ry * rz / pow(r, 5);
                    N[idx + 5] = 3. * rz * rz / pow(r, 5) - 1. / pow(r, 3);
                }
            }
        }
    }
    af::array Naf(6, n2_exp, n1_exp, n0_exp, N.data());
    Naf = reorder(Naf, 3, 2, 1, 0);
    Naf *= 1. / (4. * M_PI);
    // print("AtomisticDipoleDipoleField::N_atomistic: Naf", Naf);
    // print("Demag:Naf", Naf(0, 0, 0, af::span));
    return n2_exp == 1 ? af::fftR2C<2>(Naf) : af::fftR2C<3>(Naf);
}
} // namespace magnumaf
