#include "micro/spin_transfer_torque_zhang_li_field.hpp"
#include "constants.hpp"
#include "math.hpp"
#include "mesh.hpp"
#include "util/util.hpp"
#include <cmath>

namespace magnumaf {

SpinTransferTorqueZhangLiField::SpinTransferTorqueZhangLiField(long int j_ptr, double beta, double xi)
    : SpinTransferTorqueZhangLiField(util::pywrap::make_copy_form_py(j_ptr), beta, xi) {}

af::array field_dot_grad(af::array const& field, af::array const& m, std::array<double, 3> dx) {
    // TODO truncate on?

    // adding ghost cells like numpy 'mirror' mode
    // adding  | a b c d ... -> b | a b c d ...
    // s.t. central diff on boundary becomes zero
    auto add_ghost_cells = [](af::array const& m) {
        af::array m_ghost = af::constant(0.0, m.dims(0) + 2, m.dims(1) + 2, m.dims(2) + 2, 3, m.type());
        const auto seq_red = af::seq(1, -2);              // selecting m from m_expanded
        m_ghost(seq_red, seq_red, seq_red, af::span) = m; // setting m

        m_ghost(0, af::span, af::span, af::span) = m_ghost(2, af::span, af::span, af::span);
        m_ghost(-1, af::span, af::span, af::span) = m_ghost(-3, af::span, af::span, af::span);

        m_ghost(af::span, 0, af::span, af::span) = m_ghost(af::span, 2, af::span, af::span);
        m_ghost(af::span, -1, af::span, af::span) = m_ghost(af::span, -3, af::span, af::span);

        m_ghost(af::span, af::span, 0, af::span) = m_ghost(af::span, af::span, 2, af::span);
        m_ghost(af::span, af::span, -1, af::span) = m_ghost(af::span, af::span, -3, af::span);
        return m_ghost;
    };

    const auto m_ghost = add_ghost_cells(m);

    // derivatives of vector field along each axis (:, :, :, 3x3)
    const af::array dmdx = central_diff(m_ghost, dx[0], math::Axis::x, math::TruncateOutput::on);
    const af::array dmdy = central_diff(m_ghost, dx[1], math::Axis::y, math::TruncateOutput::on);
    const af::array dmdz = central_diff(m_ghost, dx[2], math::Axis::z, math::TruncateOutput::on);

    // Naming/extracting single vector components
    const af::array dmxdx = dmdx(af::span, af::span, af::span, 0);
    const af::array dmydx = dmdx(af::span, af::span, af::span, 1);
    const af::array dmzdx = dmdx(af::span, af::span, af::span, 2);

    const af::array dmxdy = dmdy(af::span, af::span, af::span, 0);
    const af::array dmydy = dmdy(af::span, af::span, af::span, 1);
    const af::array dmzdy = dmdy(af::span, af::span, af::span, 2);

    const af::array dmxdz = dmdz(af::span, af::span, af::span, 0);
    const af::array dmydz = dmdz(af::span, af::span, af::span, 1);
    const af::array dmzdz = dmdz(af::span, af::span, af::span, 2);

    const af::array field_x = field(af::span, af::span, af::span, 0);
    const af::array field_y = field(af::span, af::span, af::span, 1);
    const af::array field_z = field(af::span, af::span, af::span, 2);

    // applying dot grad component-wise
    const af::array field_dot_grad_x = field_x * dmxdx + field_y * dmxdy + field_z * dmxdz;
    const af::array field_dot_grad_y = field_x * dmydx + field_y * dmydy + field_z * dmydz;
    const af::array field_dot_grad_z = field_x * dmzdx + field_y * dmzdy + field_z * dmzdz;

    // rejoining components into vector field
    const af::array field_dot_grad = af::join(3, field_dot_grad_x, field_dot_grad_y, field_dot_grad_z);

    return field_dot_grad;
}

af::array field_dot_grad(af::array const& field, State const& state) {
    return field_dot_grad(field, state.m, {state.mesh.dx, state.mesh.dy, state.mesh.dz});
}

af::array SpinTransferTorqueZhangLiField::impl_H_in_Apm(const State& state) const {
    af::array Ms_field =
        state.Ms_field.isempty() ? af::constant(state.Ms, mesh::dims_s(state.mesh), state.m.type()) : state.Ms_field;
    af::array b = (beta_ * constants::mu_b) / (constants::e_abs * Ms_field * (1. + std::pow(xi_, 2.0)));
    b = af::tile(b, 1, 1, 1, 3);

    // replacing -inf and inf with zeros, handles cases in division with zero-Ms:
    // Note: arrayfire API is counterintuitive here, repalaces where values are inf (not the opposite)
    af::replace(b, !af::isInf(b), 0.0);

    af::array jgradm = field_dot_grad(j_, state);
    return (b / constants::gamma) * (math::cross4(state.m, jgradm) + xi_ * jgradm);
}
} // namespace magnumaf
