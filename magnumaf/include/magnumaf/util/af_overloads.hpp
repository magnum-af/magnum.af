#pragma once
#include "arrayfire.h"
#include "math.hpp" // vec_components
#include <stdexcept>

namespace magnumaf {

/// Overload printing af::dtype enum as readable string ( as opposed to enum/int).
inline std::ostream& operator<<(std::ostream& os, af::dtype c) {
    switch (c) {
    case f32:
        os << "f32";
        break;
    case c32:
        os << "c32";
        break;
    case f64:
        os << "f64";
        break;
    case c64:
        os << "c64";
        break;
    case b8:
        os << "b8";
        break;
    case s32:
        os << "s32";
        break;
    case u32:
        os << "u32";
        break;
    case u8:
        os << "u8";
        break;
    case s64:
        os << "s64";
        break;
    case u64:
        os << "u64";
        break;
    case s16:
        os << "s16";
        break;
    case u16:
        os << "u16";
        break;
    case f16:
        os << "f16";
        break;
    default:
        os.setstate(std::ios_base::failbit);
    }
    return os;
}

/// Prints vector components of af::array, expects size [1 1 1 3]
inline std::ostream& operator<<(std::ostream& os, const af::array& a_vector) {
    // check dimensions, expecting [1 1 1 3]
    if (a_vector.dims(0) != 1 or a_vector.dims(1) != 1 or a_vector.dims(2) != 1 or a_vector.dims(3) != 3) {
        const auto str_dims = std::to_string(a_vector.dims(0)) + ' ' + std::to_string(a_vector.dims(1)) + ' ' +
                              std::to_string(a_vector.dims(2)) + ' ' + std::to_string(a_vector.dims(3));
        const auto errormessage =
            "af::array::operator<<(): Expected af::array of dimension [1 1 1 3]. Got [" + str_dims + "].";
        throw std::runtime_error(errormessage);
    }

    switch (a_vector.type()) {
    case f64: {
        const auto [mx, my, mz] = math::vec_components<double>(a_vector);
        os << mx << '\t' << my << '\t' << mz;
    } break;
    case f32: {
        const auto [mx, my, mz] = math::vec_components<float>(a_vector);
        os << mx << '\t' << my << '\t' << mz;
    } break;
    case f16: {
        const auto [mx, my, mz] = math::vec_components<float>(a_vector.as(f32));
        os << mx << '\t' << my << '\t' << mz;
    } break;
    default: {
        const auto errormessage =
            "af::array::operator<< af::dtype::" + std::to_string(a_vector.type()) + " not supported.";
        throw std::runtime_error(errormessage);
    } break;
    }
    return os;
}

} // namespace magnumaf
