#pragma once
#include "arrayfire.h"
#include "field_terms/nonequi/nonequi_term.hpp"
#include "nonequispaced_mesh.hpp"
#include "state.hpp"
#include "util/util.hpp"
#include <functional>

namespace magnumafcpp {

class NonequiExternalField : public NonequiTerm {
  public:
    NonequiExternalField(NonequiMesh nemesh, const af::array& field) : NonequiTerm(nemesh), external_field(field) {}

    NonequiExternalField(NonequiMesh nemesh, std::function<af::array(State)> function)
        : NonequiTerm(nemesh), callback_function(function), callback_is_defined(true) {}

    ///< For wrapping only
    NonequiExternalField(NonequiMesh nemesh, long int fieldptr)
        : NonequiTerm(nemesh), external_field(util::pywrap::make_copy_form_py(fieldptr)) {}

    // Energy contribution differs by factor of 2 compared to terms linear in m
    using NonequiTerm::impl_E_in_J;
    virtual double impl_E_in_J(const State& state, const af::array& h) const override {
        return 2. * NonequiTerm::impl_E_in_J(state, h);
    };

  private:
    af::array external_field;
    std::function<af::array(State)> callback_function;
    bool callback_is_defined{false};

    virtual af::array impl_H_in_Apm(const State& state) const override {
        if (callback_is_defined) {
            return callback_function(state);
        } else {
            return external_field;
        }
    }
};

} // namespace magnumafcpp
