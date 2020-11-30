#pragma once
#include <arrayfire.h>
#include <optional>
#include <stdexcept>

namespace magnumafcpp {

/// Class for handling a variable/value which is either a double or an af::array of size [nx, ny, nz, 1]
/// overloads operators + - * / for usage.
class DoubleOrArray {
public:
  explicit DoubleOrArray(double value) : scalar(value) {}
  explicit DoubleOrArray(af::array value) : arr(value) {
      if (value.dims(3) != 1) {
          throw std::runtime_error("DoubleOrArray::DoubleOrArray(af::array): invalid input dimension, array.dims(3) != "
                                   "1. Please provide array of dimension [nx, ny, nz, 1].");
      }
  }

  // Data members, only one of the two must be value initialized
  const std::optional<double> scalar{};
  const std::optional<af::array> arr{};

  af::array operator+(const af::array& b) const {
      if (scalar) {
          // Two variants nearly equally fast:
          // return af::constant(scalar.value(), b.dims(), b.type()) + b;
          return scalar.value() + b;
      } else if (arr.value().dims() == b.dims()) {
          return arr.value() + b;
      } else if (arr.value().dims(0) == b.dims(0) and arr.value().dims(1) == b.dims(1) and
                 arr.value().dims(2) == b.dims(2) and b.dims(3) == 3) {
          return af::tile(arr.value(), 1, 1, 1, 3) + b;
      } else {
          throw std::runtime_error("DoubleOrArray::operator+: array dims do not match.");
      }
  }

  af::array operator-(const af::array& b) const {
      if (scalar) {
          return scalar.value() - b;
      } else if (arr.value().dims() == b.dims()) {
          return arr.value() - b;
      } else if (arr.value().dims(0) == b.dims(0) and arr.value().dims(1) == b.dims(1) and
                 arr.value().dims(2) == b.dims(2) and b.dims(3) == 3) {
          return af::tile(arr.value(), 1, 1, 1, 3) - b;
      } else {
          throw std::runtime_error("DoubleOrArray::operator-: array dims do not match.");
      }
  }

  af::array operator*(const af::array &b) const {
      if (scalar) {
          // return af::constant(scalar.value(), b.dims(), b.type()) * b;
          return scalar.value() * b;
      } else if (arr.value().dims() == b.dims()) {
          return arr.value() * b;
      } else if (arr.value().dims(0) == b.dims(0) and arr.value().dims(1) == b.dims(1) and
                 arr.value().dims(2) == b.dims(2) and b.dims(3) == 3) {
          return af::tile(arr.value(), 1, 1, 1, 3) * b;
      } else {
          // std::cout << "arr.value().dims()=" << arr.value().dims() << ", b.dims()=" << b.dims() <<
          // std::endl;
          throw std::runtime_error("DoubleOrArray::operator*: array dims do not match.");
      }
  }

  af::array operator/(const af::array& b) const {
      if (scalar) {
          return scalar.value() / b;
      } else if (arr.value().dims() == b.dims()) {
          return arr.value() / b;
      } else if (arr.value().dims(0) == b.dims(0) and arr.value().dims(1) == b.dims(1) and
                 arr.value().dims(2) == b.dims(2) and b.dims(3) == 3) {
          return af::tile(arr.value(), 1, 1, 1, 3) / b;
      } else {
          throw std::runtime_error("DoubleOrArray::operator/: array dims do not match.");
      }
  }

  // Getter function, returns scalar or arr as af::array
  af::array get(af::dim4 dims, af::dtype type) const {
      if (scalar) {
          return af::constant(scalar.value(), dims, type);
      } else {
          return arr.value();
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

// not commutative //inline af::array operator-(const af::array& a, const DoubleOrArray& b) { return ; }

inline af::array operator*(const af::array& a, const DoubleOrArray& b) { return b * a; }

// not commutative //inline af::array operator/(const af::array& a, const DoubleOrArray& b) { return ; }

af::array operator+(double a, const DoubleOrArray& b) = delete;
af::array operator-(double a, const DoubleOrArray& b) = delete;
af::array operator*(double a, const DoubleOrArray& b) = delete;
af::array operator/(double a, const DoubleOrArray& b) = delete;

} // namespace magnumafcpp
