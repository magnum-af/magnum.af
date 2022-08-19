#!/usr/bin/python3

import arrayfire as af
import magnumaf as maf

x, y, z = 500e-9, 125e-9, 3e-9
nx, ny, nz = 100, 25, 1

A_matr = af.constant(1.3e-11, nx, ny, nz, 1, dtype=af.Dtype.f64)

mesh = maf.Mesh(nx, ny, nz, dx=x/nx, dy=y/ny, dz=z/nz)

exc = maf.SparseExchangeField(A_matr, mesh)

# get sparse matrix
exch_matr_af = exc.get_sparse_matrix()
print("dims exch_matr_af", exch_matr_af.dims())

# convert to dense and ndarray
exch_matr = af.sparse.convert_sparse_to_dense(exch_matr_af).to_ndarray()
print("shape exch_matr", exch_matr_af.shape)
