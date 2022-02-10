#include "math.hpp"
#include <cassert>
#include <cmath>
namespace magnumaf::math {

double max_4d_abs(const af::array& a) {
    return af::max(af::max(af::max(af::max(af::abs(a), 0), 1), 2), 3).as(f64).scalar<double>();
}

af::array cross4(const af::array& a, const af::array& b) {
    af::array c = af::array(a.dims(0), a.dims(1), a.dims(2), 3, a.type());
    c(af::span, af::span, af::span, 0) = a(af::span, af::span, af::span, 1) * b(af::span, af::span, af::span, 2) -
                                         a(af::span, af::span, af::span, 2) * b(af::span, af::span, af::span, 1);
    c(af::span, af::span, af::span, 1) = a(af::span, af::span, af::span, 2) * b(af::span, af::span, af::span, 0) -
                                         a(af::span, af::span, af::span, 0) * b(af::span, af::span, af::span, 2);
    c(af::span, af::span, af::span, 2) = a(af::span, af::span, af::span, 0) * b(af::span, af::span, af::span, 1) -
                                         a(af::span, af::span, af::span, 1) * b(af::span, af::span, af::span, 0);
    return c;
}

af::array cross4shift(const af::array& a, const af::array& b) {
    af::array ashift = af::shift(a, 0, 0, 0, -1);
    af::array ashift2 = af::shift(a, 0, 0, 0, -2);
    af::array bshift = af::shift(b, 0, 0, 0, -2);
    af::array bshift2 = af::shift(b, 0, 0, 0, -1);
    return ashift * bshift - ashift2 * bshift2;
}

/// Boundary condition can by handled by providing ghost cells on expaned input array with
/// size [nx + 2, ny + 2, nz + 2, : ] and setting TruncateOutput::on, the resulting
/// af::array has reduced dims [nx, ny, nz, : ].
/// \param a input array. When incorporating BC ghost cells with dim [nx + 2, ny + 2, nz + 2, : ], set
/// TruncateOutput::on.
/// \param grid_spacing uniform grid spacing between each finite difference interval.
/// \param axis Axis enum along which to perform the finite difference, either x, y or z.
/// \param trunc_out switch wether output dims is are reduced.
///
/// Central difference is calculated along dim via
/// \f[
/// \frac{f(x_{i+1}) - f(x_{i-1})}{2h}
/// \f]
af::array central_diff(const af::array& a, double grid_spacing, Axis axis, TruncateOutput trunc_out) {
    af::array filter = af::constant(0.0, 3, 3, 3, 1, a.type());
    switch (axis) {
    case Axis::x:
        filter(0, 1, 1) = 0.5 / grid_spacing;
        filter(2, 1, 1) = -0.5 / grid_spacing;
        break;
    case Axis::y:
        filter(1, 0, 1) = 0.5 / grid_spacing;
        filter(1, 2, 1) = -0.5 / grid_spacing;
        break;
    case Axis::z:
        filter(1, 1, 0) = 0.5 / grid_spacing;
        filter(1, 1, 2) = -0.5 / grid_spacing;
        break;
    }

    auto conv = af::convolve(a, filter, AF_CONV_DEFAULT, AF_CONV_SPATIAL);

    if (trunc_out == TruncateOutput::on) {
        return conv(af::seq(1, -2), af::seq(1, -2), af::seq(1, -2), af::span); // stripping off BC ghost cells
    } else {
        return conv;
    }
}

// wrapping enum with unsigned
af::array central_diff_dimwrapped(const af::array& a, double grid_spacing, unsigned dim) {
    assert(dim < 3); // this fails
    if (dim == 0) {
        return central_diff(a, grid_spacing, Axis::x, TruncateOutput::on);
    } else if (dim == 1) {
        return central_diff(a, grid_spacing, Axis::y, TruncateOutput::on);
    } else {
        return central_diff(a, grid_spacing, Axis::z, TruncateOutput::on);
    }
}

/// 3D curl operation on vector-like af::array expanded with BC ghost cells, i.e. [nx + 2, ny + 2, nz + 2, 3]
/// Returns af::array of size [nx, ny, nz, 3]
/// \param vector_field input array [nx + 2, ny + 2, nz + 2, : ]
/// \param dx uniform grid spacing along x-axis
/// \param dy uniform grid spacing along y-axis
/// \param dz uniform grid spacing along z-axis
af::array curl_3D(const af::array& vector_field, const double dx, const double dy, const double dz,
                  TruncateOutput trunc_out) {
    assert(vector_field.dims(3) == 3);

    const auto diff_x = central_diff(vector_field, dx, Axis::x, trunc_out);
    const auto diff_y = central_diff(vector_field, dy, Axis::y, trunc_out);
    const auto diff_z = central_diff(vector_field, dz, Axis::z, trunc_out);

    auto get_x = [](const af::array& a) { return a(af::span, af::span, af::span, 0); };
    auto get_y = [](const af::array& a) { return a(af::span, af::span, af::span, 1); };
    auto get_z = [](const af::array& a) { return a(af::span, af::span, af::span, 2); };

    const auto curl_x = get_z(diff_y) - get_y(diff_z);
    const auto curl_y = get_x(diff_z) - get_z(diff_x);
    const auto curl_z = get_y(diff_x) - get_x(diff_y);

    return af::join(3, curl_x, curl_y, curl_z);
}

/// Implemented via a convolution with a FD stencil. Boundary conditions can be handled by adding ghost cells and
/// setting TruncateOutput::on
af::array laplace_3D(const af::array& a, const double dx, const double dy, const double dz, TruncateOutput trunc_out) {
    af::array filter = af::constant(0.0, 3, 3, 3, a.type());
    filter(1, 1, 1) = -2. / std::pow(dx, 2) - 2. / std::pow(dy, 2) - 2. / std::pow(dz, 2);
    filter(0, 1, 1) = 1. / std::pow(dx, 2);
    filter(2, 1, 1) = 1. / std::pow(dx, 2);
    filter(1, 0, 1) = 1. / std::pow(dy, 2);
    filter(1, 2, 1) = 1. / std::pow(dy, 2);
    filter(1, 1, 0) = 1. / std::pow(dz, 2);
    filter(1, 1, 2) = 1. / std::pow(dz, 2);
    auto conv = convolve(a, filter, AF_CONV_DEFAULT, AF_CONV_SPATIAL);
    if (trunc_out == TruncateOutput::on) {
        return conv(af::seq(1, -2), af::seq(1, -2), af::seq(1, -2), af::span); // stripping off BC ghost cells
    } else {
        return conv;
    }
}

} // namespace magnumaf::math
