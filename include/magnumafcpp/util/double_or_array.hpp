#pragma once
#include <arrayfire.h>

namespace magnumafcpp {

class DoubleOrArray {
public:
  DoubleOrArray(double value) : scalar_value(value) {}
  DoubleOrArray(af::array value) : array_value(value) {}

  // Data members
  const double scalar_value{0.};
  const af::array array_value;

  af::array operator*(const af::array &b) const {
    if (array_value.isempty()) {
      // return af::constant(scalar_value, b.dims(), b.type()) * b;
      return scalar_value * b;
    } else if (array_value.dims(3) == b.dims(3)) {
      return array_value * b;
    } else {
      return af::tile(array_value, 1, 1, 1, 3) * b;
    }
  }

  af::array operator+(const af::array &b) const {
    if (array_value.isempty()) {
      // Two variants nearly equally fast:
      // return af::constant(scalar_value, b.dims(), b.type()) + b;
      return scalar_value + b;
    } else if (array_value.dims(3) == b.dims(3)) {
      return array_value + b;
    } else {
      return af::tile(array_value, 1, 1, 1, 3) + b;
    }
  }

  // Getter function, returns scalar_value or array_value as af::array
  af::array get(af::dim4 dims, af::dtype type) const {
      if (array_value.isempty()) {
          return af::constant(scalar_value, dims, type);
      } else {
          return array_value;
      }
  }
  af::array get_as_vec(af::dim4 dims, af::dtype type) const { return af::tile(get(dims, type), 1, 1, 1, 3); }

private:
};

inline af::array operator+(const af::array &a, const DoubleOrArray &b) {
  return b + a;
}

inline af::array operator*(const af::array &a, const DoubleOrArray &b) {
  return b * a;
}

} // namespace magnumafcpp
