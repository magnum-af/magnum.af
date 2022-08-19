#include "field_terms/micro/exchange_field.hpp"
#include "math.hpp"
#include "util/util.hpp"
#include <cmath> // std::pow

namespace magnumaf {

af::array h_withedges(double A, const State& state);

// Energy calculation
// Eex=-mu0/2 integral(M . Hex) dx
// virtual parent method is overwritten as to use h_withedges
// Note: maybe this is irrelevant and can be dropped.
double ExchangeField::impl_E_in_J(const State& state, const af::array& /* h */) const {
    // TODO use h or h_widtheges?
    const auto htemp = h_withedges(this->A, state);
    if (state.Ms_field.isempty()) {
        return -constants::mu0 / 2. * state.Ms *
               af::sum(af::sum(af::sum(af::sum(htemp * state.m, 0), 1), 2), 3).scalar<double>() * state.mesh.dx *
               state.mesh.dy * state.mesh.dz;
    } else {
        return -constants::mu0 / 2. *
               af::sum(af::sum(af::sum(af::sum(state.get_Ms_field_in_vec_dims() * htemp * state.m, 0), 1), 2), 3)
                   .scalar<double>() *
               state.mesh.dx * state.mesh.dy * state.mesh.dz;
    }
}

af::array h_withedges(double A, const State& state) {
    af::array filtr = af::constant(0.0, 3, 3, 3, f64);
    // Note: skipped as this term falls out int cross product:
    filtr(1, 1, 1) =
        -2. / std::pow(state.mesh.dx, 2.) - 2. / std::pow(state.mesh.dy, 2.) - 2. / std::pow(state.mesh.dz, 2.);
    filtr(0, 1, 1) = 1 / std::pow(state.mesh.dx, 2.0);
    filtr(2, 1, 1) = 1 / std::pow(state.mesh.dx, 2.0);
    filtr(1, 0, 1) = 1 / std::pow(state.mesh.dy, 2.0);
    filtr(1, 2, 1) = 1 / std::pow(state.mesh.dy, 2.0);
    filtr(1, 1, 0) = 1 / std::pow(state.mesh.dz, 2.0);
    filtr(1, 1, 2) = 1 / std::pow(state.mesh.dz, 2.0);
    // Convolution
    af::array exch = convolve(state.m, filtr, AF_CONV_DEFAULT, AF_CONV_SPATIAL);

    // Accounting for boundary conditions by adding initial m values on the
    // boundaries by adding all 6 boundary surfaces
    exch(0, af::span, af::span, af::span) += state.m(0, af::span, af::span, af::span) / pow(state.mesh.dx, 2);
    exch(-1, af::span, af::span, af::span) += state.m(-1, af::span, af::span, af::span) / pow(state.mesh.dx, 2);

    exch(af::span, 0, af::span, af::span) += state.m(af::span, 0, af::span, af::span) / pow(state.mesh.dy, 2);
    exch(af::span, -1, af::span, af::span) += state.m(af::span, -1, af::span, af::span) / pow(state.mesh.dy, 2);

    exch(af::span, af::span, 0, af::span) += state.m(af::span, af::span, 0, af::span) / pow(state.mesh.dz, 2);
    exch(af::span, af::span, -1, af::span) += state.m(af::span, af::span, -1, af::span) / pow(state.mesh.dz, 2);

    if (state.Ms_field.isempty()) {
        return (2. * A) / (constants::mu0 * state.Ms) * exch;
    } else {
        af::array heff = (2. * A) / (constants::mu0 * state.get_Ms_field_in_vec_dims()) * exch;
        af::replace(heff, state.get_Ms_field_in_vec_dims() != 0, 0); // set all cells where Ms==0 to 0
        return heff;
    }
}

// Boundary Condition options:
enum class BC { periodic, neumann, open };

// Calculate exchange field using af::shift
af::array _shift_impl_exchange_H_in_Apm(double Aex, const State& state, BC bc) {
    const auto discrete_laplace_operator_3D_shift_impl = [](const af::array& m, double dx, double dy, double dz,
                                                            BC bc) {
        // neighbours in x+/x-/y+/y-/z+/z-:
        af::array x_plus = af::shift(m, 1) / std::pow(dx, 2.0);
        af::array x_minus = af::shift(m, -1) / std::pow(dx, 2.0);
        af::array y_plus = af::shift(m, 0, 1) / std::pow(dy, 2.0);
        af::array y_minus = af::shift(m, 0, -1) / std::pow(dy, 2.0);
        af::array z_plus = af::shift(m, 0, 0, 1) / std::pow(dz, 2.0);
        af::array z_minus = af::shift(m, 0, 0, -1) / std::pow(dz, 2.0);

        switch (bc) {
        case BC::periodic:
            // nothing to do here, this is default as generated by af::shift
            break;
        case BC::neumann:
            // TODO impl neumann BC
            throw std::runtime_error("Not yet implemented!");
            break;
        case BC::open:
            x_plus(0, af::span, af::span, af::span) = 0.0;
            x_minus(-1, af::span, af::span, af::span) = 0.0;
            y_plus(af::span, 0, af::span, af::span) = 0.0;
            y_minus(af::span, -1, af::span, af::span) = 0.0;
            z_plus(af::span, af::span, 0, af::span) = 0.0;
            z_minus(af::span, af::span, -1, af::span) = 0.0;
            break;
        }

        const af::array self = 2.0 * m / std::pow(dx, 2.0) + 2.0 * m / std::pow(dy, 2.0) + 2.0 * m / std::pow(dz, 2.0);
        af::array result = x_plus + x_minus + y_plus + y_minus + z_plus + z_minus - self;
        return result;
    };

    const af::array dlaplace =
        discrete_laplace_operator_3D_shift_impl(state.m, state.mesh.dx, state.mesh.dy, state.mesh.dz, bc);
    if (state.Ms_field.isempty()) {
        af::array H_in_Apm = 2.0 * Aex / (constants::mu0 * state.Ms) * dlaplace;
        return H_in_Apm;
    } else {
        af::array H_in_Apm = 2.0 * Aex / (constants::mu0 * state.get_Ms_field_in_vec_dims()) * dlaplace;
        af::replace(H_in_Apm, state.get_Ms_field_in_vec_dims() != 0, 0); // set all cells where Ms==0 to 0
        return H_in_Apm;
    }
}

// Terms proportional to m drop out in the cross product of the LLG and thus is
// neglected as arrayfire is extremely slow with indexing operations NOTE: This
// yields no longer the physical exchange field but optimizes the caluclation
af::array _convolve_old_impl_exchange_H_in_Apm(double A, const State& state) {
    af::array filtr = af::constant(0.0, 3, 3, 3, f64);
    // Note: skippable as this term falls out int cross product:
    filtr(1, 1, 1) =
        -2. / std::pow(state.mesh.dx, 2.) - 2. / std::pow(state.mesh.dy, 2.) - 2. / std::pow(state.mesh.dz, 2.);
    filtr(0, 1, 1) = 1 / pow(state.mesh.dx, 2);
    filtr(2, 1, 1) = 1 / pow(state.mesh.dx, 2);
    filtr(1, 0, 1) = 1 / pow(state.mesh.dy, 2);
    filtr(1, 2, 1) = 1 / pow(state.mesh.dy, 2);
    filtr(1, 1, 0) = 1 / pow(state.mesh.dz, 2);
    filtr(1, 1, 2) = 1 / pow(state.mesh.dz, 2);
    af::array exch = convolve(state.m, filtr, AF_CONV_DEFAULT, AF_CONV_SPATIAL);
    if (state.Ms_field.isempty()) {
        return (2. * A) / (constants::mu0 * state.Ms) * exch;
    } else {
        af::array heff = (2. * A) / (constants::mu0 * state.get_Ms_field_in_vec_dims()) * exch;
        af::replace(heff, state.get_Ms_field_in_vec_dims() != 0, 0); // set all cells where Ms==0 to 0
        return heff;
    }
}

af::array _convolve_lap_impl_exchange_H_in_Apm(double A, const State& state) {
    const auto lapl = math::laplace_3D(state.m, state.mesh.dx, state.mesh.dy, state.mesh.dz, math::TruncateOutput::off);
    if (state.Ms_field.isempty()) {
        return (2. * A) / (constants::mu0 * state.Ms) * lapl;
    } else {
        af::array heff = (2. * A) / (constants::mu0 * state.get_Ms_field_in_vec_dims()) * lapl;
        af::replace(heff, state.get_Ms_field_in_vec_dims() != 0, 0); // set all cells where Ms==0 to 0
        return heff;
    }
}

af::array ExchangeField::impl_H_in_Apm(const State& state) const {
    if constexpr (false) { // TODO enable as only option
        return _shift_impl_exchange_H_in_Apm(this->A, state, BC::open);
    } else {
        if constexpr (true) {
            return _convolve_old_impl_exchange_H_in_Apm(this->A, state);
        } else {
            return _convolve_lap_impl_exchange_H_in_Apm(this->A, state);
        }
    }
}
} // namespace magnumaf
