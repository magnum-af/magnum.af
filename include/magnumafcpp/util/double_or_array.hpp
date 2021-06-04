#pragma once
#include "util/variant_helper.hpp"
#include <arrayfire.h>
#include <stdexcept>
#include <variant>

namespace magnumafcpp::util {

/// Class for handling a variable/value which is either a double or an af::array of size [nx, ny, nz, 1]
/// overloads operators + - * / for usage.
class DoubleOrArray {
  public:
    explicit DoubleOrArray(double value) : variant_(value) {}

    // Note: af::array ctor is not explicit, so any ctor with int/double/... can be mapped to this ctor
    // This can cause nasty implicit-conversion bugs, e.g. 'double 5' calls ctor for af::array(5) and not double!
    explicit DoubleOrArray(af::array value) : variant_(value) {
        if (value.dims(3) != 1) {
            throw std::runtime_error(
                "DoubleOrArray::DoubleOrArray(af::array): invalid input dimension, array.dims(3) != "
                "1. Please provide array of dimension [nx, ny, nz, 1].");
        }
    }

    af::array operator+(const af::array& b) const {
        struct Add {
            const af::array& b;
            af::array operator()(double d) { return d + b; }
            af::array operator()(const af::array& a) {
                if (a.dims() == b.dims()) {
                    return a + b;
                } else if (a.dims(0) == b.dims(0) and a.dims(1) == b.dims(1) and a.dims(2) == b.dims(2) and
                           b.dims(3) == 3) {
                    return af::tile(a, 1, 1, 1, 3) + b;
                } else {
                    throw std::runtime_error("DoubleOrArray::operator+: array dims do not match.");
                }
            }
        };
        return std::visit(Add{b}, variant_);
    }

    af::array operator-(const af::array& b) const {
        const auto visitor =
            make_visitor([&b](double d) { return d - b; },
                         [&b](const af::array& a) {
                             if (a.dims() == b.dims()) {
                                 return a - b;
                             } else if (a.dims(0) == b.dims(0) and a.dims(1) == b.dims(1) and a.dims(2) == b.dims(2) and
                                        b.dims(3) == 3) {
                                 return af::tile(a, 1, 1, 1, 3) - b;
                             } else {
                                 throw std::runtime_error("DoubleOrArray::operator-: array dims do not match.");
                             }
                         });
        return std::visit(visitor, variant_);
    }

    af::array operator*(const af::array& b) const {
        const auto visitor = overloaded{
            [&b](double d) { return d * b; },
            [&b](const af::array& a) {
                if (a.dims() == b.dims()) {
                    return a * b;
                } else if (a.dims(0) == b.dims(0) and a.dims(1) == b.dims(1) and a.dims(2) == b.dims(2) and
                           b.dims(3) == 3) {
                    return af::tile(a, 1, 1, 1, 3) * b;
                } else {
                    throw std::runtime_error("DoubleOrArray::operator*: array dims do not match.");
                }
            },
        };
        return std::visit(visitor, variant_);
    }

    af::array operator/(const af::array& b) const {
        const auto visitor = overloaded{
            [&b](double d) { return d / b; },
            [&b](const af::array& a) {
                const auto divide = [&b](const af::array& a) {
                    if (a.dims() == b.dims()) {
                        return a / b;
                    } else if (a.dims(0) == b.dims(0) and a.dims(1) == b.dims(1) and a.dims(2) == b.dims(2) and
                               b.dims(3) == 3) {
                        return af::tile(a, 1, 1, 1, 3) / b;
                    } else {
                        throw std::runtime_error("DoubleOrArray::operator/: array dims do not match.");
                    }
                };
                auto result = divide(a);
                af::replace(result, b != 0, 0); // replacing potential divs by null (i.e. NaN) with zeros.
                                                // Optional optimization: cache whether a has zero vals
                return result;
            },
        };
        return std::visit(visitor, variant_);
    }

    // Getter function, returns stored value as af::array
    af::array operator()(af::dim4 dims, af::dtype type) const {
        const auto visitor = overloaded{
            [&](double d) { return af::constant(d, dims, type); },
            [&](const af::array& a) {
                if (a.dims() == dims) {
                    return a.as(type);
                } else if (a.dims(0) == dims[0] and a.dims(1) == dims[1] and a.dims(2) == dims[2] and dims[3] == 3) {
                    return af::tile(a, 1, 1, 1, 3).as(type);
                } else {
                    throw std::runtime_error("DoubleOrArray::operator(): array dims do not match.");
                }
            },
        };
        return std::visit(visitor, variant_);
    }

  private:
    std::variant<double, af::array> variant_; // Note: be careful when adding variant types which can be implicitly
                                              // converted to af::array (s.a. int/float)

    // Preventing erroneous operators
    // double/int/.. b deduces af::array(b), which creates an [b,1,1,1] array type f32 type
    af::array operator+(double b) const = delete;
    af::array operator-(double b) const = delete;
    af::array operator*(double b) const = delete;
    af::array operator/(double b) const = delete;
};

inline af::array get_as_vec(const DoubleOrArray& a, af::dim4 dims, af::dtype type) {
    return af::tile(a(dims, type), 1, 1, 1, 3);
}

inline af::array operator+(const af::array& a, const DoubleOrArray& b) { return b + a; }

// not commutative
inline af::array operator-(const af::array& a, const DoubleOrArray& b) {
    const auto btemp = b(a.dims(), a.type());
    return a - btemp;
}

inline af::array operator*(const af::array& a, const DoubleOrArray& b) { return b * a; }

// not commutative
inline af::array operator/(const af::array& a, const DoubleOrArray& b) {
    const auto btemp = b(a.dims(), a.type());
    af::array result = a / btemp;
    af::replace(result, btemp != 0, 0);
    return result;
}

af::array operator+(double a, const DoubleOrArray& b) = delete;
af::array operator-(double a, const DoubleOrArray& b) = delete;
af::array operator*(double a, const DoubleOrArray& b) = delete;
af::array operator/(double a, const DoubleOrArray& b) = delete;

} // namespace magnumafcpp::util
