#pragma once
#include "LLGTerm.hpp"
#include "../state.hpp"
#include "arrayfire.h"
#include <functional>

namespace magnumaf{


class ExternalField : public LLGTerm {
  public:
    ExternalField(af::array zee_in);///< Constant Zeeman field.
    ExternalField(af::array (*callback_func_in)(State state));///< Callback function for e.g. time dependent external field
    ExternalField(std::function<af::array(State)>);
    ExternalField(long int zee_in_addr);///< For wrapping only

    af::array h(const State& state);//Field contribution
    long int h_ptr(const State& state);// For wrapping

    float E(const State& state);//Energy contribution
    float E(const State& state, const af::array& h);///< Calculating the micromagnetic energy for a already calculated h field

    long int get_m_addr();// For wrapping only
    void set_homogeneous_field(const float x, const float y, const float z);///< Setting homogeneous zeeman field with x, y, z components of static Zeeman field.

    float get_cpu_time(){return 0;}// use or remove

  private:
    af::array calc_heff(const State& state);
    af::array zee_field;
    af::array (*callback_func)(State state);
    const bool callback{false};
    std::function<af::array(State)> lamda_callback;
    const bool is_lamda{false};
    //float af_time{0.};
};


//for wrapping:
//https://stackoverflow.com/questions/8800838/how-to-pass-a-function-pointer-to-an-external-program-in-cython
}// namespace magnumaf
