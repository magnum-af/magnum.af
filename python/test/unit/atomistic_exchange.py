import unittest
import arrayfire as af
import numpy as np
import magnumaf
import math

class AtomisticExchangeFieldTest(unittest.TestCase):
  # arbitrary parameters:
  p  = 1e-12
  J  = 1e-12
  dx = 2.715e-10
  #TODO# add init test#self.assertEqual(self.J_atom, J_atom)

  def test_atomistic_exchange_2_1_1_z_z(self):
    mesh=magnumaf.Mesh(2, 1, 1, self.dx, self.dx, self.dx)
    p =self.p
    J_atom =self.J
    m=af.constant(0.0, 2, 1, 1, 3, dtype=af.Dtype.f64)
    m[0, 0, 0, 0] = 0
    m[0, 0, 0, 1] = 0
    m[0, 0, 0, 2] = 1

    m[1, 0, 0, 0] = 0
    m[1, 0, 0, 1] = 0
    m[1, 0, 0, 2] = 1

    state=magnumaf.State(mesh, Ms = p, m = m)
    atom_exch=magnumaf.AtomisticExchangeField(J_atom)
    llg=magnumaf.LLGIntegrator(alpha = 0, terms = [atom_exch])

    self.assertAlmostEqual(llg.E(state), -J_atom)

    af_heff = llg.h(state)
    np_heff = af_heff.__array__()

    self.assertAlmostEqual(np_heff[0, 0, 0, 0], 0 )
    self.assertAlmostEqual(np_heff[0, 0, 0, 1], 0 )
    self.assertAlmostEqual(np_heff[0, 0, 0, 2], J_atom/magnumaf.Constants.mu0/p )
    self.assertAlmostEqual(np_heff[1, 0, 0, 0], 0 )
    self.assertAlmostEqual(np_heff[1, 0, 0, 1], 0 )
    self.assertAlmostEqual(np_heff[1, 0, 0, 2], J_atom/magnumaf.Constants.mu0/p )

  def test_atomistic_exchange_2_1_1_z_x(self):
    mesh=magnumaf.Mesh(2, 1, 1, self.dx, self.dx, self.dx)
    p =self.p
    J_atom =self.J
    m=af.constant(0.0, 2, 1, 1, 3, dtype=af.Dtype.f64)
    m[0, 0, 0, 0] = 0
    m[0, 0, 0, 1] = 0
    m[0, 0, 0, 2] = 1

    m[1, 0, 0, 0] = 1
    m[1, 0, 0, 1] = 0
    m[1, 0, 0, 2] = 0

    state=magnumaf.State(mesh, Ms = p, m = m)
    atom_exch=magnumaf.AtomisticExchangeField(J_atom)
    llg=magnumaf.LLGIntegrator(alpha = 0, terms = [atom_exch])

    self.assertAlmostEqual(llg.E(state), 0)

    af_heff = llg.h(state)
    np_heff = af_heff.__array__()

    self.assertAlmostEqual(np_heff[0, 0, 0, 0], J_atom/magnumaf.Constants.mu0/p )
    self.assertAlmostEqual(np_heff[0, 0, 0, 1], 0 )
    self.assertAlmostEqual(np_heff[0, 0, 0, 2], 0 )
    self.assertAlmostEqual(np_heff[1, 0, 0, 0], 0 )
    self.assertAlmostEqual(np_heff[1, 0, 0, 1], 0 )
    self.assertAlmostEqual(np_heff[1, 0, 0, 2], J_atom/magnumaf.Constants.mu0/p )

  def test_atomistic_exchange_2_1_1_z_minusz(self):
    mesh=magnumaf.Mesh(2, 1, 1, self.dx, self.dx, self.dx)
    p =self.p
    J_atom =self.J
    m=af.constant(0.0, 2, 1, 1, 3, dtype=af.Dtype.f64)
    m[0, 0, 0, 0] = 0
    m[0, 0, 0, 1] = 0
    m[0, 0, 0, 2] = 1

    m[1, 0, 0, 0] = 0
    m[1, 0, 0, 1] = 0
    m[1, 0, 0, 2] =-1

    state=magnumaf.State(mesh, Ms = p, m = m)
    atom_exch=magnumaf.AtomisticExchangeField(J_atom)
    llg=magnumaf.LLGIntegrator(alpha = 0, terms = [atom_exch])

    self.assertAlmostEqual(llg.E(state), J_atom)

    af_heff = llg.h(state)
    np_heff = af_heff.__array__()

    self.assertAlmostEqual(np_heff[0, 0, 0, 0], 0 )
    self.assertAlmostEqual(np_heff[0, 0, 0, 1], 0 )
    self.assertAlmostEqual(np_heff[0, 0, 0, 2], -J_atom/magnumaf.Constants.mu0/p )
    self.assertAlmostEqual(np_heff[1, 0, 0, 0], 0 )
    self.assertAlmostEqual(np_heff[1, 0, 0, 1], 0 )
    self.assertAlmostEqual(np_heff[1, 0, 0, 2], J_atom/magnumaf.Constants.mu0/p )

  def test_atomistic_exchange_1_2_1_z_z(self):
    mesh=magnumaf.Mesh(1, 2, 1, self.dx, self.dx, self.dx)
    p =self.p
    J_atom =self.J
    m=af.constant(0.0, 1, 2, 1, 3, dtype=af.Dtype.f64)
    m[0, 0, 0, 0] = 0
    m[0, 0, 0, 1] = 0
    m[0, 0, 0, 2] = 1

    m[0, 1, 0, 0] = 0
    m[0, 1, 0, 1] = 0
    m[0, 1, 0, 2] = 1

    state=magnumaf.State(mesh, Ms = p, m = m)
    atom_exch=magnumaf.AtomisticExchangeField(J_atom)
    llg=magnumaf.LLGIntegrator(alpha = 0, terms = [atom_exch])

    self.assertAlmostEqual(llg.E(state), -J_atom)

    af_heff = llg.h(state)
    np_heff = af_heff.__array__()

    self.assertAlmostEqual(np_heff[0, 0, 0, 0], 0 )
    self.assertAlmostEqual(np_heff[0, 0, 0, 1], 0 )
    self.assertAlmostEqual(np_heff[0, 0, 0, 2], J_atom/magnumaf.Constants.mu0/p )
    self.assertAlmostEqual(np_heff[0, 1, 0, 0], 0 )
    self.assertAlmostEqual(np_heff[0, 1, 0, 1], 0 )
    self.assertAlmostEqual(np_heff[0, 1, 0, 2], J_atom/magnumaf.Constants.mu0/p )

if __name__ == '__main__':
  unittest.main()
