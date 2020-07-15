#pragma once
#include "arrayfire.h"

namespace magnumafcpp {

std::ostream& operator<<(std::ostream& os, af::dtype c) {
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

} // namespace magnumafcpp
