#pragma once
#include "arrayfire.h"
#include <memory>

namespace magnumafcpp::util {

/// Wrapper for a std::unique_ptr<T[], > with custom deleter to read af::array raw data.
/// Typename T must match array type. If a.isempty() ptr is set to nullptr.
/// Note: dtor handels case ptr==nullptr by not calling get_deleter() when get()==nullptr.
template <typename T> class HostPtrAccessor {
  public:
    HostPtrAccessor(const af::array& a) : ptr(a.isempty() ? nullptr : a.host<T>(), &af::freeHost) {}
    const T& operator[](size_t index) const { return ptr[index]; }
    explicit operator bool() const noexcept { return ptr ? true : false; }

  private:
    std::unique_ptr<T[], decltype(&af::freeHost)> ptr{};
};

template <typename T> class HostPtrManipulator {
  public:
    HostPtrManipulator(af::array& a) : ptr(a.isempty() ? nullptr : a.host<T>(), &af::freeHost) {}
    T& operator[](size_t index) const { return ptr[index]; }
    explicit operator bool() const noexcept { return ptr ? true : false; }

  private:
    std::unique_ptr<T[], decltype(&af::freeHost)> ptr{};
};

/// RAII wrapper for raw data read access of af::array
/// Note nullptr handling is not straight-forward, therefore prefer using unique_ptr
template <typename T> class HostRawPtrAccessor {
  public:
    HostRawPtrAccessor(const af::array& a) : ptr(a.host<T>()) {}
    ~HostRawPtrAccessor() { af::freeHost(ptr); }
    const T& operator[](size_t index) const { return ptr[index]; }

  private:
    const T* const ptr;
};

/// RAII wrapper for raw data read/write access of af::array
template <typename T> class HostRawPtrManipulator {
  public:
    HostRawPtrManipulator(af::array& a) : ptr(a.host<T>()) {}
    ~HostRawPtrManipulator() { af::freeHost(ptr); }
    T& operator[](size_t index) const { return ptr[index]; }

  private:
    T* const ptr;
};

} // namespace magnumafcpp::util
