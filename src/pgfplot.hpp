#ifndef PGFPLOT_H
#define PGFPLOT_H
#include <iostream>
#include <fstream>
#include "arrayfire.h"
#include "mesh.hpp"
#include "func.hpp"

void pgfplot_mz(af::array m, Mesh mesh, std::string outputfile, int mz_slice);
void pgfplot_nz_quiver(af::array m, Mesh mesh, std::string outputfile, int mz_slice, int nmod);

#endif
