#pragma once
#include "arrayfire.h"

namespace magnumafcpp::util {

/// RAII wrapper for raw data read access of af::array
template <typename T> class HostPtrAccessor {
  public:
    HostPtrAccessor(const af::array& a) : ptr(a.host<T>()) {}
    ~HostPtrAccessor() { af::freeHost(ptr); }
    const T& operator[](size_t index) const { return ptr[index]; }

  private:
    const T* const ptr;
};

/// RAII wrapper for raw data read/write access of af::array
template <typename T> class HostPtrManipulator {
  public:
    HostPtrManipulator(af::array& a) : ptr(a.host<T>()) {}
    ~HostPtrManipulator() { af::freeHost(ptr); }
    T& operator[](size_t index) const { return ptr[index]; }

  private:
    T* const ptr;
};

} // namespace magnumafcpp::util
