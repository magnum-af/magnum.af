#pragma once
#include <arrayfire.h>
#include <stdexcept>

namespace magnumafcpp {

class DoubleOrArray {
public:
  DoubleOrArray(double value) : scalar_value(value) {}
  DoubleOrArray(af::array value) : array_value(value) {
      if (value.dims(3) != 1) {
          throw std::runtime_error("DoubleOrArray::DoubleOrArray(af::array): invalid input dimension, array.dims(3) != "
                                   "1. Please provide array of dimension [nx, ny, nz, 1].");
      }
  }

  // Data members
  const double scalar_value{0.};
  const af::array array_value;

  af::array operator+(const af::array& b) const {
      if (array_value.isempty()) {
          // Two variants nearly equally fast:
          // return af::constant(scalar_value, b.dims(), b.type()) + b;
          return scalar_value + b;
      } else if (array_value.dims() == b.dims()) {
          return array_value + b;
      } else if (array_value.dims(0) == b.dims(0) and array_value.dims(1) == b.dims(1) and
                 array_value.dims(2) == b.dims(2) and b.dims(3) == 3) {
          return af::tile(array_value, 1, 1, 1, 3) + b;
      } else {
          throw std::runtime_error("DoubleOrArray::operator+: array dims do not match.");
      }
  }

  af::array operator-(const af::array& b) const {
      if (array_value.isempty()) {
          return scalar_value - b;
      } else if (array_value.dims() == b.dims()) {
          return array_value - b;
      } else if (array_value.dims(0) == b.dims(0) and array_value.dims(1) == b.dims(1) and
                 array_value.dims(2) == b.dims(2) and b.dims(3) == 3) {
          return af::tile(array_value, 1, 1, 1, 3) - b;
      } else {
          throw std::runtime_error("DoubleOrArray::operator-: array dims do not match.");
      }
  }

  af::array operator*(const af::array &b) const {
      if (array_value.isempty()) {
          // return af::constant(scalar_value, b.dims(), b.type()) * b;
          return scalar_value * b;
      } else if (array_value.dims() == b.dims()) {
          return array_value * b;
      } else if (array_value.dims(0) == b.dims(0) and array_value.dims(1) == b.dims(1) and
                 array_value.dims(2) == b.dims(2) and b.dims(3) == 3) {
          return af::tile(array_value, 1, 1, 1, 3) * b;
      } else {
          // std::cout << "array_value.dims()=" << array_value.dims() << ", b.dims()=" << b.dims() << std::endl;
          throw std::runtime_error("DoubleOrArray::operator*: array dims do not match.");
      }
  }

  af::array operator/(const af::array& b) const {
      if (array_value.isempty()) {
          return scalar_value / b;
      } else if (array_value.dims() == b.dims()) {
          return array_value / b;
      } else if (array_value.dims(0) == b.dims(0) and array_value.dims(1) == b.dims(1) and
                 array_value.dims(2) == b.dims(2) and b.dims(3) == 3) {
          return af::tile(array_value, 1, 1, 1, 3) / b;
      } else {
          throw std::runtime_error("DoubleOrArray::operator/: array dims do not match.");
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
  // Preventing erroneous operators
  // double/int/.. b deduces af::array(b), which creates an [b,1,1,1] array type f32 type
  af::array operator+(double b) const = delete;
  af::array operator-(double b) const = delete;
  af::array operator*(double b) const = delete;
  af::array operator/(double b) const = delete;
};

inline af::array operator+(const af::array& a, const DoubleOrArray& b) { return b + a; }

inline af::array operator-(const af::array& a, const DoubleOrArray& b) { return b - a; }

inline af::array operator*(const af::array& a, const DoubleOrArray& b) { return b * a; }

inline af::array operator/(const af::array& a, const DoubleOrArray& b) { return b / a; }

af::array operator+(double a, const DoubleOrArray& b) = delete;
af::array operator-(double a, const DoubleOrArray& b) = delete;
af::array operator*(double a, const DoubleOrArray& b) = delete;
af::array operator/(double a, const DoubleOrArray& b) = delete;

} // namespace magnumafcpp
