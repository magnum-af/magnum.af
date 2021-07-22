#pragma once
#include <utility>
namespace magnumafcpp::util {

// helper constant for type-matching visitor (#3) (https://en.cppreference.com/w/cpp/utility/variant/visit)
template <class> inline constexpr bool always_false_v = false;
// overload pattern (#4)
template <class... Ts> struct overloaded : Ts... { using Ts::operator()...; };
// explicit deduction guide (not needed as of C++20, g++10.2 still needs it)
template <class... Ts> overloaded(Ts...) -> overloaded<Ts...>;

/// make visitor from lambdas (c++20):
#if (__cplusplus > 201703L)
template <class... Ts> auto make_visitor(Ts&&... ts) { return overloaded<Ts...>(std::forward<Ts>(ts)...); }
#else
// make_visitor (c++11) (adapted from https://bitbashing.io/std-visit.html):
template <class... Fs> struct mk_vis_overload;
template <class F0, class... Frest> struct mk_vis_overload<F0, Frest...> : F0, mk_vis_overload<Frest...> {
    mk_vis_overload(F0 f0, Frest... rest) : F0(f0), mk_vis_overload<Frest...>(rest...) {}
    using F0::operator();
    using mk_vis_overload<Frest...>::operator();
};
template <class F0> struct mk_vis_overload<F0> : F0 {
    mk_vis_overload(F0 f0) : F0(f0) {}
    using F0::operator();
};
template <class... Fs> auto make_visitor(Fs&&... fs) { return mk_vis_overload<Fs...>(std::forward<Fs>(fs)...); }
#endif

} // namespace magnumafcpp::util
