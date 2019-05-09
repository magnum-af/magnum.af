#ifndef ZEE_H
#define ZEE_H
#include <functional>
#include "arrayfire.h"
#include "LLGTerm.hpp"
#include "../state.hpp"
#include "../func.hpp"

class ExternalField : public LLGTerm {
  public:
    ///< Constant Zeeman field.
    ExternalField(af::array zee_in);
    ExternalField(long int zee_in_addr);
    long int get_m_addr();
    ///< Setting x,y,z components of static Zeeman field.
    void set_xyz(const State&, const double x, const double y, const double z);

    ///< Time dependent Zeeman field.
    ExternalField(af::array (*callback_func_in)(State state));
    ExternalField(std::function<af::array(State)>);
    af::array (*callback_func)(State state);
    std::function<af::array(State)> lamda_callback;
    bool callback{false};
    bool is_lamda{false};

    //Field contribution
    af::array h(const State& state);
    //Energy contribution
    double E(const State& state);
    double E(const State& state, const af::array& h);///< Calculating the micromagnetic energy for a already calculated h field
    double get_cpu_time(){return cpu_time;}//TODO remove or use
    double cpu_time{0.};
    af::array zee_field;
    af::timer timer;
};


#endif

//for wrapping: 
//https://stackoverflow.com/questions/8800838/how-to-pass-a-function-pointer-to-an-external-program-in-cython
