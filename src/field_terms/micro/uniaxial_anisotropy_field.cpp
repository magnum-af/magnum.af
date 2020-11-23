// Ref Diss Abert p 15 sec 2.3 eq 2.3.28
// Hu(r)=2 Ku1 /(mu0 Ms) eu ( eu . m)
// With higher order (not implemented): Hu(r)=2 Ku1 /(mu0 Ms) eu ( eu . m) ( + 4
// Ku2 /(mu0 Ms) eu ( eu . m)^3
#include "micro/uniaxial_anisotropy_field.hpp"
#include "util/func.hpp"
#include "util/misc.hpp" // Warning()
#include "util.hpp"

namespace magnumafcpp {

UniaxialAnisotropyField::UniaxialAnisotropyField(double Ku1, std::array<double, 3> Ku1_axis)
    : Ku1(Ku1), Ku1_axis(normalize_vector(Ku1_axis)) {}

UniaxialAnisotropyField::UniaxialAnisotropyField(af::array Ku1_field, std::array<double, 3> Ku1_axis)
    : Ku1_field(Ku1_field.dims(3) == 1 ? af::tile(Ku1_field, 1, 1, 1, 3) : Ku1_field),
      Ku1_axis(normalize_vector(Ku1_axis)) {
    if (Ku1_field.dims(3) == 3) {
        printf("%s UniaxialAnisotropyField: You are using legacy dimension "
               "[nx, ny, nz, 3] for Ku1, please now use scalar field "
               "dimensions [nx, ny, nz, 1].\n",
               Warning());
    }
}

UniaxialAnisotropyField::UniaxialAnisotropyField(af::array Ku1_field, af::array Ku1_axis_field)
    : Ku1_field(Ku1_field.dims(3) == 1 ? af::tile(Ku1_field, 1, 1, 1, 3) : Ku1_field), Ku1_axis_field(Ku1_axis_field) {
    if (Ku1_field.dims(3) == 3) {
        printf("%s UniaxialAnisotropyField: You are using legacy dimension "
               "[nx, ny, nz, 3] for Ku1, please now use scalar field "
               "dimensions [nx, ny, nz, 1].\n",
               Warning());
    }
}

UniaxialAnisotropyField::UniaxialAnisotropyField(double Ku1, long int Ku1_axis_field_ptr)
    : Ku1(Ku1), Ku1_axis_field(*(new af::array(*((void**)Ku1_axis_field_ptr)))) {}

// For wrapping only
UniaxialAnisotropyField::UniaxialAnisotropyField(long int Ku1_field_ptr, long int Ku1_axis_field_ptr)
    : Ku1_field((*(new af::array(*((void**)Ku1_field_ptr)))).dims(3) == 1
                    ? af::tile(*(new af::array(*((void**)Ku1_field_ptr))), 1, 1, 1, 3)
                    : *(new af::array(*((void**)Ku1_field_ptr)))),
      Ku1_axis_field(*(new af::array(*((void**)Ku1_axis_field_ptr)))) {
    if ((*(new af::array(*((void**)Ku1_field_ptr)))).dims(3) == 3) {
        printf("%s UniaxialAnisotropyField You are using legacy dimension [nx, "
               "ny, nz, 3] for Ku1, please now use scalar field dimensions "
               "[nx, ny, nz, 1].\n",
               Warning());
    }
}

// For wrapping only
UniaxialAnisotropyField::UniaxialAnisotropyField(double Ku1, double Ku1_axis_0, double Ku1_axis_1, double Ku1_axis_2)
    : Ku1(Ku1), Ku1_axis(normalize_vector(std::array<double, 3>{Ku1_axis_0, Ku1_axis_1, Ku1_axis_2})) {}

// For wrapping only
UniaxialAnisotropyField::UniaxialAnisotropyField(long int Ku1_field_ptr, double Ku1_axis_0, double Ku1_axis_1,
                                                 double Ku1_axis_2)
    : Ku1_field((*(new af::array(*((void**)Ku1_field_ptr)))).dims(3) == 1
                    ? af::tile(*(new af::array(*((void**)Ku1_field_ptr))), 1, 1, 1, 3)
                    : *(new af::array(*((void**)Ku1_field_ptr)))),
      Ku1_axis(normalize_vector(std::array<double, 3>{Ku1_axis_0, Ku1_axis_1, Ku1_axis_2})) {
    if ((*(new af::array(*((void**)Ku1_field_ptr)))).dims(3) == 3) {
        printf("%s UniaxialAnisotropyField: You are using legacy dimension "
               "[nx, ny, nz, 3] for Ku1, please now use scalar field "
               "dimensions [nx, ny, nz, 1].\n",
               Warning());
    }
}

af::array UniaxialAnisotropyField::h(const State& state) {
    af::timer timer_anisotropy = af::timer::start();

    // switch Ku1_axis and Ku1_axis_field
    af::array eu; // Array containing normal vectors
    if (Ku1_axis_field.isempty()) {
        eu = af::array(dims_vector(state.mesh), f64);
        eu(af::span, af::span, af::span, 0) = Ku1_axis[0];
        eu(af::span, af::span, af::span, 1) = Ku1_axis[1];
        eu(af::span, af::span, af::span, 2) = Ku1_axis[2];
    } else {
        eu = Ku1_axis_field;
    }

    af::array anisotropy = eu * state.m;
    anisotropy = af::sum(anisotropy, 3);
    anisotropy = af::tile(anisotropy, 1, 1, 1, 3);

    if (state.afsync)
        af::sync();
    computation_time_heff += af::timer::stop(timer_anisotropy);
    if (state.Ms_field.isempty() && Ku1_field.isempty()) {
        return 2. * Ku1 / (constants::mu0 * state.Ms) * (eu * anisotropy);
    } else if (!state.Ms_field.isempty() && Ku1_field.isempty()) {
        af::array result = 2. * Ku1 / (constants::mu0 * state.Ms_field) * (eu * anisotropy);
        af::replace(result, !af::iszero(state.Ms_field),
                    0); // Replacing all resulting NaN with 0
        return result;
    } else if (state.Ms_field.isempty() && !Ku1_field.isempty()) {
        return 2. * Ku1_field / (constants::mu0 * state.Ms) * (eu * anisotropy);
    } else {
        af::array result = 2. * Ku1_field / (constants::mu0 * state.Ms_field) * (eu * anisotropy);
        af::replace(result, !af::iszero(state.Ms_field),
                    0); // Replacing all resulting NaN with 0
        return result;
    }
}

double UniaxialAnisotropyField::get_ku1_axis(int i) { return Ku1_axis[i]; }

// void UniaxialAnisotropyField::set_Ku1_field(long int aptr){
//    void **a = (void **)aptr;
//    micro_Ku1_field = *( new af::array( *a ));
//}

long int UniaxialAnisotropyField::get_Ku1_field() {
    af::array* a = new af::array(Ku1_field);
    return (long int)a->get();
}
} // namespace magnumafcpp
