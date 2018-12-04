#ifndef ZEE_H
#define ZEE_H
#include <functional>
#include "arrayfire.h"
#include "LLGTerm.hpp"
#include "../state.hpp"
#include "../func.hpp"

class Zee : public LLGTerm {
  public:
    //static Zeeman field
    Zee(af::array zee_in);
    Zee(long int zee_in_addr);

    //time dependent Zeeman field
    Zee(af::array (*callback_func_in)(State state));
    Zee(std::function<af::array(State)>);
    af::array (*callback_func)(State state);
    std::function<af::array(State)> lamda_callback;
    bool callback{false};
    bool is_lamda{false};

    //Field contribution
    af::array h(const State& state);
    //Energy contribution
    double E(const State& state);
    double E(const State& state, const af::array& h);///< Calculating the micromagnetic energy for a already calculated h field
    //CPU time
    double get_cpu_time(){return cpu_time;}

    double rate;//[T/s]
    double hzee_max;//[T]
    af::array zee_field;
    double cpu_time{0.};
    af::timer timer;
};


#endif

//for wrapping: 
//https://stackoverflow.com/questions/8800838/how-to-pass-a-function-pointer-to-an-external-program-in-cython
