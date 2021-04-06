#pragma once
#include "arrayfire.h"
#include <array>
/// Common math functions for af::array container.
namespace magnumafcpp::math{

/// Calculate mean along first three dimensions of af::array.
template <typename T> std::array<T, 3> mean_3d(const af::array& a) {
    const auto mean = af::mean(af::mean(af::mean(a, 0), 1), 2);
    const T mx = mean(0, 0, 0, 0).scalar<T>();
    const T my = mean(0, 0, 0, 1).scalar<T>();
    const T mz = mean(0, 0, 0, 2).scalar<T>();
    return {mx, my, mz};
}

/// Absolute value of maximum of all values in array
double max_4d_abs(const af::array& a);

/// Cross product for vector fields of format [nx, ny, nz, 3]
af::array cross4(const af::array& a, const af::array& b);

/// Cross product for vector fields of format [nx, ny, nz, 3] using af::shift
af::array cross4shift(const af::array& a, const af::array& b); // Note: Slightly slower than cross4

/// Axis along which to perform operation s.a. central_diff.
enum class Axis { x, y, z };

/// Truncate ghost boundary cells introduced due to boundary conditions.
/// This reduces the size of the output by 2 in x,y and z.
enum class TruncateOutput { on, off };

/// Central finite difference, computed via convolution.
af::array central_diff(const af::array& a, double grid_spacing, Axis, TruncateOutput);

/// Curl vector operation.
af::array curl_3D(const af::array& vector_field, const double dx, const double dy, const double dz, TruncateOutput);

/// Laplace operation.
af::array laplace_3D(const af::array& a, const double dx, const double dy, const double dz, TruncateOutput);

} // namespace magnumafcpp::math
