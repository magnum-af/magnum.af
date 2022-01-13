#pragma once
#include "util/util.hpp"
#include "util/variant_helper.hpp"
#include <array>
#include <arrayfire.h>
#include <variant>

namespace magnumaf::util {

///
/// Wrapper for a unit vector either of type std::array<double, 3> or af::array dims [nx, ny, nz, 3]
///
struct UnitVectorOrArray {
    explicit UnitVectorOrArray(std::array<double, 3> scalar_vector) : variant(util::normalize_vector(scalar_vector)) {}

    // Note: af::array ctor is not explicit, so any further variants s.a. with int/double/... can be mapped to this ctor
    explicit UnitVectorOrArray(af::array array_vector) : variant(util::normalize(array_vector)) {
        if (array_vector.dims(3) != 3) {
            throw std::runtime_error("VectorOrArray::VectorOrArray(af::array): invalid input dimension, array.dims(3) "
                                     "!= 3. Please provide array of dimension [nx, ny, nz, 3].");
        }
    }

    std::variant<std::array<double, 3>, af::array> variant;

    /// Getter function, returns array or tiles vector to af::array with dimension dims
    /// \param dims expects af::dim4(nx, ny, nz) to tile the unit vector, optional last (4th) dimension is neglected.
    /// \param type sets af::dtype of returned af::array
    af::array get_as_array(af::dim4 dims, af::dtype type) const {
        const auto visitor = [&](auto&& array) -> af::array {
            using T = std::decay_t<decltype(array)>;
            if constexpr (std::is_same_v<T, std::array<double, 3>>) {
                const af::array vector_as_array = af::array(1, 1, 1, 3, array.data());
                return af::tile(vector_as_array.as(type), dims.dims[0], dims.dims[1], dims.dims[2], 1);
            } else if constexpr (std::is_same_v<T, af::array>) {
                return array.as(type);
            } else {
                static_assert(always_false_v<T>, "non-exhaustive visitor!");
            }
        };
        return std::visit(visitor, variant);
    }
};
} // namespace magnumaf::util
