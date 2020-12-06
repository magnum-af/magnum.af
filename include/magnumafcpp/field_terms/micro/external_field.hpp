#pragma once
#include "state.hpp"
#include "field_term.hpp"
#include "arrayfire.h"
#include <functional>

namespace magnumafcpp {

class ExternalField : public FieldTerm {
  public:
    ExternalField(af::array zee_in);                           ///< Constant Zeeman field.
    ExternalField(af::array (*callback_func_in)(State state)); ///< Callback function for e.g. time dependent external
                                                               ///< field
    ExternalField(std::function<af::array(State)>);
    ExternalField(long int zee_in_addr); ///< For wrapping only

    virtual af::array h(const State& state) const override; // Field contribution

    using FieldTerm::E;
    virtual double E(const State& state, const af::array& h) const override;

    void set_homogeneous_field(const double x, const double y,
                               const double z); ///< Setting homogeneous zeeman field with x, y, z
                                                ///< components of static Zeeman field.

  private:
    af::array zee_field;
    af::array (*callback_func)(State state);
    const bool callback{false};
    std::function<af::array(State)> lamda_callback;
    const bool is_lamda{false};
    // double af_time{0.};
};

// for wrapping:
// https://stackoverflow.com/questions/8800838/how-to-pass-a-function-pointer-to-an-external-program-in-cython
} // namespace magnumafcpp
