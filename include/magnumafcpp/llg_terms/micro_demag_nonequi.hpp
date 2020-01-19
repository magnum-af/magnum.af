#pragma once
#include "LLGTerm.hpp"
#include "../state.hpp"
#include "arrayfire.h"

namespace magnumafcpp
{

namespace newell_nonequi
{
double Nxx(const double x, const double y, const double z, const double dx, const double dy, const double dz, const double dX, const double dY, const double dZ);
double Nxy(const double x, const double y, const double z, const double dx, const double dy, const double dz, const double dX, const double dY, const double dZ);
double nonequi_index_distance(const std::vector<double> spacings, const unsigned i, const unsigned j, const bool verbose = true);
} // namespace newell_nonequi

class NonEquiDemagField : public LLGTerm
{
public:
    //Field contribution
    af::array h(const State &state);
    //Energy contribution
    double E(const State &state);
    double E(const State &state, const af::array &h); ///< Calculating the micromagnetic energy for a already calculated h field
    //CPU time
    double get_cpu_time() { return cpu_time; }

    NonEquiDemagField(NonequispacedMesh nonequimesh, bool verbose = true, bool caching = false, unsigned nthreads = 0);

    af::array Nfft; //!< Array storing the Fourier transfrom of the demag tensor.

    double cpu_time{0.};
    af::timer timer_demagsolve;

    //For wrapping
    void print_Nfft();

    //std::vector<double> z_spacing;
private:
    const unsigned nthreads;
    af::array calculate_N(int n0_exp, int n1_exp, int n2_exp, double dx, double dy, const std::vector<double> z_spacing);
};
} // namespace magnumafcpp
