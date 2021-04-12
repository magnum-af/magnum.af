#pragma once
namespace magnumafcpp::util {

template <class> inline constexpr bool always_false_v = false;
// helper type for the visitor #4
template <class... Ts> struct overloaded : Ts... { using Ts::operator()...; };
// explicit deduction guide (not needed as of C++20)
template <class... Ts> overloaded(Ts...) -> overloaded<Ts...>;

} // namespace magnumafcpp::util
