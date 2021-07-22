#include "field_terms/atom/atomistic_external_field.hpp"

namespace magnumafcpp {

double AtomisticExternalField::impl_E_in_J(const State& state, const af::array& h) const {
    if (state.Ms_field.isempty()) {
        return -constants::mu0 * state.Ms * sum(sum(sum(sum(h * state.m, 0), 1), 2), 3).scalar<double>();
    } else {
        return -constants::mu0 * sum(sum(sum(sum(state.Ms_field * h * state.m, 0), 1), 2), 3).scalar<double>();
    }
}

} // namespace magnumafcpp