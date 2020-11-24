#pragma once
#include "arrayfire.h"

namespace magnumafcpp::util {

template <typename T> struct HostPtrHandle {
    HostPtrHandle(const af::array& a) : ptr(a.host<T>()) {}
    ~HostPtrHandle() { af::freeHost(ptr); }
    const T* ptr;
};
} // namespace magnumafcpp::util
