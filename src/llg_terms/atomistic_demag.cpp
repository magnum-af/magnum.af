#include "atomistic_demag.hpp"
#include "../func.hpp"

namespace magnumafcpp {

af::array N_atomistic(int n0_exp, int n1_exp, int n2_exp, double dx, double dy,
                      double dz);

AtomisticDipoleDipoleField::AtomisticDipoleDipoleField(Mesh mesh) {
    Nfft = N_atomistic(mesh.n0_exp, mesh.n1_exp, mesh.n2_exp, mesh.dx, mesh.dy,
                       mesh.dz);
    // h   =af::array (mesh.n0_exp    , mesh.n1_exp, mesh.n2_exp, 3, f64);
}

af::array AtomisticDipoleDipoleField::h(const State& state) {
    timer_demagsolve = af::timer::start();
    // FFT with zero-padding of the m field
    af::array mfft;
    if (state.mesh.n2_exp == 1) {
        mfft = af::fftR2C<2>(state.m,
                             af::dim4(state.mesh.n0_exp, state.mesh.n1_exp));
    } else {
        mfft = af::fftR2C<3>(
            state.m,
            af::dim4(state.mesh.n0_exp, state.mesh.n1_exp, state.mesh.n2_exp));
    }

    af::array hfft = af::array(state.mesh.n0_exp / 2 + 1, state.mesh.n1_exp,
                               state.mesh.n2_exp, 3, c64);
    // Pointwise product
    hfft(af::span, af::span, af::span, 0) =
        Nfft(af::span, af::span, af::span, 0) *
            mfft(af::span, af::span, af::span, 0) +
        Nfft(af::span, af::span, af::span, 1) *
            mfft(af::span, af::span, af::span, 1) +
        Nfft(af::span, af::span, af::span, 2) *
            mfft(af::span, af::span, af::span, 2);
    hfft(af::span, af::span, af::span, 1) =
        Nfft(af::span, af::span, af::span, 1) *
            mfft(af::span, af::span, af::span, 0) +
        Nfft(af::span, af::span, af::span, 3) *
            mfft(af::span, af::span, af::span, 1) +
        Nfft(af::span, af::span, af::span, 4) *
            mfft(af::span, af::span, af::span, 2);
    hfft(af::span, af::span, af::span, 2) =
        Nfft(af::span, af::span, af::span, 2) *
            mfft(af::span, af::span, af::span, 0) +
        Nfft(af::span, af::span, af::span, 4) *
            mfft(af::span, af::span, af::span, 1) +
        Nfft(af::span, af::span, af::span, 5) *
            mfft(af::span, af::span, af::span, 2);

    // IFFT reversing padding
    af::array h_field;
    if (state.mesh.n2_exp == 1) {
        h_field = af::fftC2R<2>(hfft);
        // af::print("h_dip", state.Ms * h_field(seq(0, state.mesh.n0_exp/2-1),
        // seq(0, state.mesh.n1_exp/2-1)));//TODO hack
        if (state.afsync)
            af::sync();
        cpu_time += af::timer::stop(timer_demagsolve);
        return state.Ms *
               h_field(af::seq(0, state.mesh.n0_exp / 2 - 1),
                       af::seq(0, state.mesh.n1_exp / 2 -
                                      1)); // TODO consider p density, then we
                                           // have to multip at m before fft
    } else {
        h_field = af::fftC2R<3>(hfft);
        if (state.afsync)
            af::sync();
        cpu_time += af::timer::stop(timer_demagsolve);
        return state.Ms * h_field(af::seq(0, state.mesh.n0_exp / 2 - 1),
                                  af::seq(0, state.mesh.n1_exp / 2 - 1),
                                  af::seq(0, state.mesh.n2_exp / 2 - 1),
                                  af::span);
    }
}

af::array N_atomistic(int n0_exp, int n1_exp, int n2_exp, double dx, double dy,
                      double dz) {
    double* N = NULL;
    N = new double[n0_exp * n1_exp * n2_exp * 6];
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
    af::array Naf(6, n2_exp, n1_exp, n0_exp, N);
    Naf = reorder(Naf, 3, 2, 1, 0);
    Naf *= 1. / (4. * M_PI);
    // print("AtomisticDipoleDipoleField::N_atomistic: Naf", Naf);
    // print("Demag:Naf", Naf(0, 0, 0, af::span));
    delete[] N;
    N = NULL;
    if (n2_exp == 1) {
        Naf = af::fftR2C<2>(Naf);
    } else {
        Naf = af::fftR2C<3>(Naf);
    }
    return Naf;
}
} // namespace magnumafcpp
