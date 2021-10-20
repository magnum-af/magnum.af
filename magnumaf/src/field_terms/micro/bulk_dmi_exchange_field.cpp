#include "field_terms/micro/bulk_dmi_exchange_field.hpp"
#include "math.hpp"
#include "util/color_string.hpp"
#include <array>
#include <iostream>

namespace magnumaf {

BulkDMIExchangeField::BulkDMIExchangeField(double D, double A) : D_bulk_(D), A_(A) {}

af::array get_m_with_BC_ghost_cells(double D, double A, const af::array& m, double dx, double dy, double dz) {
    const auto x = 0; // for readability of last dim
    const auto y = 1; // for readability of last dim
    const auto z = 2; // for readability of last dim
    // would also work, could cause unexpected behaviour: enum Vec{x, y, z};
    const auto seq_red = af::seq(1, -2); // selecting m from m_expanded

    // Note: af::array() would cause random values on cells not overwritten.
    af::array m_with_ghost_cells = af::constant(0.0, m.dims(0) + 2, m.dims(1) + 2, m.dims(2) + 2, 3, m.type());

    m_with_ghost_cells(seq_red, seq_red, seq_red, af::span) = m; // setting m

    // Boundary direction along axis, +x would be plus; -x minus
    // Internally used to choose indices and correct sign
    enum class Boundary { plus, minus };

    auto set_x = [D, A, &m_with_ghost_cells, m, seq_red](double dxyz, Boundary boundary) {
        const int i = boundary == Boundary::plus ? -1 : 0;       // index at boundary
        const int i_off_1 = boundary == Boundary::plus ? -2 : 1; // index at boundary one cell inward

        const auto coeffs = [D_bulk_ = D, A_ = A, dxyz]() { return 2.0 * dxyz * D_bulk_ / (2.0 * A_); };
        const auto mx_at_i_off_1 = m(i_off_1, af::span, af::span, x);
        const auto my_at_i_off_1 = m(i_off_1, af::span, af::span, y);
        const auto mz_at_i_off_1 = m(i_off_1, af::span, af::span, z);

        const double sign = boundary == Boundary::plus ? 1.0 : -1.0;
        const auto minus_mz_at_i = sign * -1. * m(i, af::span, af::span, z);
        const auto my_at_i = sign * m(i, af::span, af::span, y);

        m_with_ghost_cells(i, seq_red, seq_red, x) = mx_at_i_off_1;                            // dmx/dx
        m_with_ghost_cells(i, seq_red, seq_red, y) = my_at_i_off_1 + coeffs() * minus_mz_at_i; // dmy/dx
        m_with_ghost_cells(i, seq_red, seq_red, z) = mz_at_i_off_1 + coeffs() * my_at_i;       // dmz/dx
    };

    auto set_y = [D, A, &m_with_ghost_cells, m, seq_red](double dxyz, Boundary boundary) {
        const int i = boundary == Boundary::plus ? -1 : 0;       // index at boundary
        const int i_off_1 = boundary == Boundary::plus ? -2 : 1; // index at boundary one cell inward

        auto coeffs = [D_bulk_ = D, A_ = A, dxyz]() { return 2.0 * dxyz * D_bulk_ / (2.0 * A_); };

        const auto mx_at_i_off_1 = m(af::span, i_off_1, af::span, x);
        const auto my_at_i_off_1 = m(af::span, i_off_1, af::span, y);
        const auto mz_at_i_off_1 = m(af::span, i_off_1, af::span, z);

        const double sign = boundary == Boundary::plus ? 1.0 : -1.0;
        const auto mz_at_i = sign * m(af::span, i, af::span, z);
        const auto minus_mx_at_i = sign * -1. * m(af::span, i, af::span, x);

        m_with_ghost_cells(seq_red, i, seq_red, x) = mx_at_i_off_1 + coeffs() * mz_at_i;       // dmx/dx
        m_with_ghost_cells(seq_red, i, seq_red, y) = my_at_i_off_1;                            // dmx/dy
        m_with_ghost_cells(seq_red, i, seq_red, z) = mz_at_i_off_1 + coeffs() * minus_mx_at_i; // dmx/dz
    };

    auto set_z = [D, A, &m_with_ghost_cells, m, seq_red](double dxyz, Boundary boundary) {
        const int i = boundary == Boundary::plus ? -1 : 0;       // index at boundary
        const int i_off_1 = boundary == Boundary::plus ? -2 : 1; // index at boundary one cell inward

        auto coeffs = [D_bulk_ = D, A_ = A, dxyz]() { return 2.0 * dxyz * D_bulk_ / (2.0 * A_); };

        const auto mx_at_i_off_1 = m(af::span, af::span, i_off_1, x);
        const auto my_at_i_off_1 = m(af::span, af::span, i_off_1, y);
        const auto mz_at_i_off_1 = m(af::span, af::span, i_off_1, z);

        const double sign = boundary == Boundary::plus ? 1.0 : -1.0;
        const auto minus_my_at_i = sign * -1. * m(af::span, af::span, i, y);
        const auto mx_at_i = sign * m(af::span, af::span, i, x);

        m_with_ghost_cells(seq_red, seq_red, i, x) = mx_at_i_off_1 + coeffs() * minus_my_at_i;
        m_with_ghost_cells(seq_red, seq_red, i, y) = my_at_i_off_1 + coeffs() * mx_at_i;
        m_with_ghost_cells(seq_red, seq_red, i, z) = mz_at_i_off_1;
    };

    set_x(dx, Boundary::minus); // at -x, i.e. i=0
    set_x(dx, Boundary::plus);  // at +x, i.e. i=nx-1

    set_y(dy, Boundary::minus); // at -y
    set_y(dy, Boundary::plus);  // at +y

    set_z(dz, Boundary::minus); // at -z
    set_z(dz, Boundary::plus);  // at +z

    return m_with_ghost_cells;
}

af::array BulkDMIExchangeField::impl_H_in_Apm(const State& state) const {
    if (state.m.dims(0) < 3 or state.m.dims(1) < 3 or state.m.dims(2) < 3) {
        std::cout << color_string::warning()
                  << "BulkDMIExchangeField: number of cells in x, y or z < 3. Special 2D case is not yet handled, "
                     "still applying BC ghost cells."
                  << std::endl;
    }

    const auto m_with_ghost_cells =
        get_m_with_BC_ghost_cells(D_bulk_, A_, state.m, state.mesh.dx, state.mesh.dy, state.mesh.dz);
    const auto buld_dmi =
        2.0 * D_bulk_ *
        math::curl_3D(m_with_ghost_cells, state.mesh.dx, state.mesh.dy, state.mesh.dz, math::TruncateOutput::on);

    const auto exchange =
        -2.0 * A_ *
        math::laplace_3D(m_with_ghost_cells, state.mesh.dx, state.mesh.dy, state.mesh.dz, math::TruncateOutput::on);
    const auto result = buld_dmi + exchange;

    if (state.Ms_field.isempty()) {
        return result / (-constants::mu0 * state.Ms);
    } else {
        return result / (-constants::mu0 * state.get_Ms_field_in_vec_dims());
    }
}

} // namespace magnumaf
