#pragma once
#include "arrayfire.h"
#include "field_terms/nonequi/nonequi_term.hpp"
#include "state.hpp"

namespace magnumaf {

class NonequiDemagField : public NonequiTerm {
  public:
    NonequiDemagField(const NonequiMesh& nonequimesh, bool verbose = true, bool caching = false, unsigned nthreads = 0);

  private:
    af::array Nfft; //!< Array storing the Fourier transfrom of the demag tensor.
    virtual af::array impl_H_in_Apm(const State& state) const override; // Field contribution
};

// expose for testing:
namespace newell_nonequi {
double Nxx(const double x, const double y, const double z, const double dx, const double dy, const double dz,
           const double dX, const double dY, const double dZ);
double Nxy(const double x, const double y, const double z, const double dx, const double dy, const double dz,
           const double dX, const double dY, const double dZ);
double nonequi_index_distance(const std::vector<double>& spacings, const unsigned i, const unsigned j,
                              const bool verbose = true);

} // namespace newell_nonequi
} // namespace magnumaf
