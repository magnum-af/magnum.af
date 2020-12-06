#pragma once
#include "func.hpp"
#include "util.hpp"
#include <array>
#include <arrayfire.h>
#include <optional>

namespace magnumafcpp {

///
/// Wrapper for a unit vector either of type std::array<double, 3> or af::array dims [nx, ny, nz, 3]
///
class UnitVectorOrArray {
  public:
    const std::optional<std::array<double, 3>> scalar_vector;
    const std::optional<af::array> array_vector; // dims [nx, ny, nz, 3]

    UnitVectorOrArray(std::array<double, 3> scalar_vector) : scalar_vector(normalize_vector(scalar_vector)) {}

    UnitVectorOrArray(af::array array_vector) : array_vector(normalize_handle_zero_vectors(array_vector)) {
        if (array_vector.dims(3) != 3) {
            throw std::runtime_error("VectorOrArray::VectorOrArray(af::array): invalid input dimension, array.dims(3) "
                                     "!= 3. Please provide array of dimension [nx, ny, nz, 3].");
        }
    }

    /// Getter function, returns array or tiles vector to af::array
    /// arugments dims and type are only used if (scalar_vector)
    /// \param dims expects [nx, ny, nz, X] to tile the unit vector to, last dimension is neglected.
    /// Only used if (scalar_vector)
    /// \param type sets af::dtype of returning af::array
    af::array get_as_array(af::dim4 dims, af::dtype type) const {
        if (scalar_vector) {
            af::array vector_as_array = af::array(1, 1, 1, 3, scalar_vector.value().data());
            return af::tile(vector_as_array.as(type), dims.dims[0], dims.dims[1], dims.dims[2], 1);
        } else {
            if (dims != array_vector.value().dims()) {
            }
            return array_vector.value().as(type);
        }
    }

  private:
};
} // namespace magnumafcpp
