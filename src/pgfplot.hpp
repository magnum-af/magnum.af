#pragma once
#include "mesh.hpp"
#include "arrayfire.h"

namespace magnumaf{


void pgfplot_mz(af::array m, Mesh mesh, std::string outputfile, int mz_slice);
void pgfplot_mz(af::array m, Mesh mesh, std::string outputfile, int mz_slice, int nmod);
void pgfplot_mi(af::array m, Mesh mesh, std::string outputfile, int i, int mz_slice, int nmod);
void pgfplot_nz_quiver(af::array m, Mesh mesh, std::string outputfile, int mz_slice, int nmod);

}// namespace magnumaf
