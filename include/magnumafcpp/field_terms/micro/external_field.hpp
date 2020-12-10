#pragma once
#include "arrayfire.h"
#include "field_term.hpp"
#include "state.hpp"
#include <functional>

namespace magnumafcpp {

class ExternalField : public FieldTerm {
  public:
    explicit ExternalField(af::array zee_in);                           ///< Constant Zeeman field.
    explicit ExternalField(af::array (*callback_func_in)(State state)); ///< Callback function for e.g. time dependent
                                                                        ///< external field
    explicit ExternalField(std::function<af::array(State)>);
    explicit ExternalField(long int zee_in_addr); ///< For wrapping only

    void set_homogeneous_field(const double x, const double y,
                               const double z); ///< Setting homogeneous zeeman field with x, y, z
                                                ///< components of static Zeeman field.

  private:
    af::array zee_field;
    af::array (*callback_func)(State state);
    bool callback{false};
    std::function<af::array(State)> lamda_callback;
    bool is_lamda{false};
    // double af_time{0.};

    virtual af::array impl_H_in_Apm(const State& state) const override; // Field contribution

    virtual double impl_E_in_J(const State& state, const af::array& h) const override;
};

// for wrapping:
// https://stackoverflow.com/questions/8800838/how-to-pass-a-function-pointer-to-an-external-program-in-cython
} // namespace magnumafcpp
