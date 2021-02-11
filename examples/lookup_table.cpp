// Demo script on using a lookup table for setting material parameters and creating the full array on the fly

#include "arrayfire.h"
#include "util/af_overloads.hpp"
#include <array>
#include <iostream>
#include <vector>

int main() {
    using namespace magnumafcpp; // Note: gets << operator overload into scope

    {
        std::vector<double> _values = {0.5, 1.5, 2.5, 3.5};
        af::array values(_values.size(), _values.data());
        af::print("values", values);
        std::cout << "values.type()=" << values.type() << std::endl;

        std::vector<uint8_t> _index = {2, 1, 0, 3, 4}; // Note: 4 i out of bound an defaults to last index
        af::array index(_index.size(), _index.data());
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
        af::array values(_values.size(), _values.data());
        af::print("values", values);
        std::cout << "values.type()=" << values.type() << std::endl;

        const auto dims = af::dim4(2, 4, 6, 1);
        af::array index = af::constant(0, dims, u8);
        index(af::span, af::span, af::seq(0, (index.dims(2) - 1) / 2)) = 1;
        index(af::span, af::seq(0, (index.dims(1) - 1) / 2), af::seq(0, (index.dims(2) - 1) / 2)) = 2;
        index(af::seq(0, (index.dims(0) - 1) / 2), af::seq(0, (index.dims(1) - 1) / 2),
              af::seq(0, (index.dims(2) - 1) / 2)) = 3;
        af::print("index", index);
        std::cout << "index.type()=" << index.type() << std::endl;

        auto make_lookup = [](const af::array& values, const af::array& index) {
            return af::moddims(af::lookup(values, af::flat(index)), index.dims());
        };
        af::array lookup = make_lookup(values, index);

        af::print("lookup", lookup);
        std::cout << "lookup.type()=" << lookup.type() << std::endl;

        auto max_abs_diff = [](const af::array& a, const af::array& b) {
            return af::max(af::max(af::max(af::max(af::abs(a - b), 0), 1), 2), 3);
        };
        std::cout << "max abs diff should be zero if values equal index: max_abs_diff = "
                  << max_abs_diff(lookup, index).scalar<double>() << std::endl;
    }

    return 0;
}
