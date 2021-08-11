// Demo script on using a lookup table for setting material parameters and creating the full array on the fly
// Notes: we should exclusively use af::array for values of lookup (at least internally) to preserve dynamic type
// For convenience, we can provide templated ctors which convert std containers to af::array
// Fixing interface to e.g. std::vector<double> would fix values to af::dtype::f64

#include "arrayfire.h"
#include "mesh.hpp"
#include "util/af_overloads.hpp"
#include <array>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <vector>

inline af::array do_multidim_lookup(const af::array& values, const af::array& index) {
    return af::moddims(af::lookup(values, af::flat(index)), index.dims());
}

template <typename T> af::array make_afarray_from_container(T const& container) {
    if (container.size() == 0) {
        throw std::runtime_error("do_lookup: values.size() is zero! Must be >= 1.");
    }
    return af::array(container.size(), container.data());
}

// Lookup on optional region.
// If region not set, uses backup_dims for return size
inline af::array apply_lookup_on_optional(af::array const& values, std::optional<af::array> const& opt_region,
                                          af::dim4 const& backup_dims) {
    if (opt_region) {
        return do_multidim_lookup(values, opt_region.value());
    } else {
        return do_multidim_lookup(values, af::constant(0, backup_dims, u8));
    }
}

// Idea: unique m, region (= index), A = {A0, A1, A2, ...}
// State has m, Ms and region, field_terms only have vecs with values
// TODO: how to handle m, Ms ODR for zero values?
struct RegionState {
    magnumaf::Mesh mesh;
    af::array m;
    af::array Ms;
    std::optional<af::array> regions;
};

// lookup wrapper for RegionState.
// We use free function with wrapper to keeping logic more seperated from data structure
inline af::array apply_lookup(af::array const& values, RegionState const& region_state) {
    return apply_lookup_on_optional(values, region_state.regions, magnumaf::mesh::dims_s(region_state.mesh));
}

// convenience wrapper
template <typename T> auto apply_lookup_from_container(T const& values, RegionState const& region_state) {
    return apply_lookup(make_afarray_from_container(values), region_state);
}

// exemplary interaction which holds af::array values internally, not std::vector or such
struct Interaction {
    af::array values{};
    explicit Interaction(af::array values) : values(std::move(values)) {}
    // template <typename T> Interaction(T const& values) : Interaction(make_afarray_from_container(values)) {}
    // convenience ctor for passing std containers:
    // NOTE: could/should be constrained to std::is_floating_point, i.e. float, double as std::is_integral types like
    // int, uint, char not wanted in most cases
    template <typename T>
    explicit Interaction(std::vector<T> const& values) : Interaction(make_afarray_from_container(values)) {}
    // convenience ctor for passing single values like float, double:
    // Note: using vector overload instead of af::array i.o.t. deduce af::array type
    template <typename T> explicit Interaction(T const& value) : Interaction(std::vector<T>{value}) {}
};

inline af::array apply_lookup(Interaction const& inter, RegionState const& region_state) {
    return apply_lookup(inter.values, region_state);
}

int main() {
    using namespace magnumaf; // Note: gets << operator overload into scope

    {
        std::vector<double> _values = {0.5, 1.5, 2.5, 3.5};
        af::array values = make_afarray_from_container(_values);
        af::print("values", values);
        std::cout << "values.type()=" << values.type() << std::endl;

        std::vector<uint8_t> _index = {2, 1, 0, 3, 4}; // Note: 4 i out of bound an defaults to last index
        af::array index(make_afarray_from_container(_index));
        af::print("index", index);
        std::cout << "index.type()=" << index.type() << std::endl;

        af::print("index", index);
        af::array lookup = af::lookup(values, index);
        af::print("lookup", lookup);
        std::cout << "lookup.type()=" << lookup.type() << std::endl;
    }

    // use arbitrarily sized index array for lookup:
    {
        constexpr std::array<double, 4> _values = {0, 1, 2, 3};
        // Note: for dynamic input use: std::vector<double> _values = {...};
        af::array values = make_afarray_from_container(_values);
        // af::array values(_values.size(), _values.data());
        af::print("values", values);
        std::cout << "values.type()=" << values.type() << std::endl;

        const auto dims = af::dim4(2, 4, 6, 1);
        af::array index = af::constant(0, dims, u8);
        index(af::span, af::span, af::seq(0, (index.dims(2) - 1) / 2)) = 1;
        index(af::span, af::seq(0, (index.dims(1) - 1) / 2), af::seq(0, (index.dims(2) - 1) / 2)) = 2;
        index(af::seq(0, (index.dims(0) - 1) / 2), af::seq(0, (index.dims(1) - 1) / 2),
              af::seq(0, (index.dims(2) - 1) / 2)) = 3;
        af::print("index 2", index);
        std::cout << "index.type()=" << index.type() << std::endl;

        af::array lookup = do_multidim_lookup(values, index);

        af::print("lookup 2", lookup);
        std::cout << "lookup.type()=" << lookup.type() << std::endl;

        auto max_abs_diff = [](const af::array& a, const af::array& b) {
            return af::max(af::max(af::max(af::max(af::abs(a - b), 0), 1), 2), 3);
        };
        std::cout << "max abs diff should be zero if values equal index: max_abs_diff = "
                  << max_abs_diff(lookup, index).scalar<double>() << std::endl;
    }

    {
        // defining mesh and region
        Mesh mesh{3, 1, 1, 1e-9, 2e-9, 3e-9};
        auto dims_scal = mesh::dims_s(mesh);
        af::array region = af::constant(0, dims_scal, af::dtype::u8);
        for (uint8_t i = 0; i < mesh.nx; ++i) {
            region(i, af::span, af::span, af::span) = i;
        }
        af::print("region", region);

        // using an af::array region for lookup
        auto print_array_and_type0 = [i = 0](auto const& values, af::array const& region) mutable {
            const auto result = do_multidim_lookup(make_afarray_from_container(values), region);
            af::print(("lookup" + std::to_string(i)).c_str(), result);
            std::cout << result.type() << std::endl;
            i++;
        };
        print_array_and_type0(std::vector{0.5, 1.5, 2.5}, region);
        print_array_and_type0(std::vector{0, 1, -2}, region);   // implicit <int>
        print_array_and_type0(std::vector{0U, 1U, 2U}, region); // <unsigned>

        // using RegionState for lookup
        auto dims_vec = mesh::dims_s(mesh);
        af::array m = af::constant(0, dims_vec, f64);
        m(af::span, af::span, af::span, 0) = 1;
        af::array Ms = af::constant(0, dims_scal, af::dtype::f64);
        const auto region_state = RegionState{mesh, m, Ms, region};

        const auto print_array_and_type = [&rs = std::as_const(region_state)](auto const& values) {
            const auto result = apply_lookup_from_container(values, rs);
            af::print("", result);
            std::cout << result.type() << std::endl;
        };

        print_array_and_type(std::vector{0.5, 1.5, 2.5}); // <double>
        print_array_and_type(std::vector{-1, 0, 1, 2});   // <int>
        print_array_and_type(std::vector{0U, 1U, 2U});    // unsigned

        const auto noregion_state = RegionState{mesh, m, Ms, {}};
        try {
            af::print("noregion_state", apply_lookup_from_container(std::vector<double>{}, noregion_state));
        } catch (std::exception const& e) {
            std::cerr << e.what() << std::endl;
        }

        const auto make_inter = [&rs = std::as_const(region_state)](auto vals) {
            Interaction inter{vals};
            const auto result = apply_lookup(inter, rs);
            af::print("make_inter:", result);
            std::cout << result.type() << std::endl;
        };

        make_inter(af::iota(af::dim4(mesh.nx), af::dim4(1), f64));
        make_inter(af::iota(af::dim4(mesh.nx), af::dim4(1), u32));
        make_inter(af::iota(af::dim4(mesh.nx), af::dim4(1), f32));
        make_inter(std::vector{0.5, 1.5, 4.5});
        make_inter(std::vector{1, 2, 3});
        make_inter(2.5);
        make_inter(1);
    }

    return 0;
}
