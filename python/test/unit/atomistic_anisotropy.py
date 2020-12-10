import unittest
import arrayfire as af
import numpy as np
import magnumaf
import math

class AtomisticAnisotropyTest(unittest.TestCase):
  # arbitrary parameters:
  p  = 1e-20
  k  = 1e-20
  dx = 2.715e-10
  #TODO init test
  #TODO#Ku1_atom_axis = [0., 3., 0.]
    #TODO#self.assertEqual((0., 1., 0.), Ku1_atom_axis)

  def test_atomistic_anisotropy_2_1_1_z_z(self):
    mesh=magnumaf.Mesh(2, 1, 1, self.dx, self.dx, self.dx)
    p=self.p
    Ku1_atom=self.k
    m=af.constant(0.0, 2, 1, 1, 3, dtype=af.Dtype.f64)
    m[0, 0, 0, 0] = 0
    m[0, 0, 0, 1] = 0
    m[0, 0, 0, 2] = 1

    m[1, 0, 0, 0] = 0
    m[1, 0, 0, 1] = 0
    m[1, 0, 0, 2] = 1

    state=magnumaf.State(mesh, Ms = p, m = m)
    atom_ani=magnumaf.AtomisticUniaxialAnisotropyField(Ku1_atom, [0, 0, 1])

    self.assertAlmostEqual(atom_ani.Energy_in_J(state), -2*Ku1_atom)

    af_heff = atom_ani.h(state)
    np_heff = af_heff.__array__()

    self.assertAlmostEqual(np_heff[0, 0, 0, 0], 0 )
    self.assertAlmostEqual(np_heff[0, 0, 0, 1], 0 )
    self.assertAlmostEqual(np_heff[0, 0, 0, 2], 2*Ku1_atom/magnumaf.Constants.mu0/p )
    self.assertAlmostEqual(np_heff[1, 0, 0, 0], 0 )
    self.assertAlmostEqual(np_heff[1, 0, 0, 1], 0 )
    self.assertAlmostEqual(np_heff[1, 0, 0, 2], 2*Ku1_atom/magnumaf.Constants.mu0/p )

  def test_atomistic_anisotropy_2_1_1_z_x(self):
    mesh=magnumaf.Mesh(2, 1, 1, self.dx, self.dx, self.dx)
    p =self.p
    Ku1_atom =self.k
    m=af.constant(0.0, 2, 1, 1, 3, dtype=af.Dtype.f64)
    m[0, 0, 0, 0] = 0
    m[0, 0, 0, 1] = 0
    m[0, 0, 0, 2] = 1

    m[1, 0, 0, 0] = 1
    m[1, 0, 0, 1] = 0
    m[1, 0, 0, 2] = 0

    state=magnumaf.State(mesh, Ms = p, m = m)
    atom_ani=magnumaf.AtomisticUniaxialAnisotropyField(Ku1_atom, [0, 0, 1])

    self.assertAlmostEqual(atom_ani.Energy_in_J(state), -Ku1_atom)

    af_heff = atom_ani.h(state)
    np_heff = af_heff.__array__()

    self.assertAlmostEqual(np_heff[0, 0, 0, 0], 0 )
    self.assertAlmostEqual(np_heff[0, 0, 0, 1], 0 )
    self.assertAlmostEqual(np_heff[0, 0, 0, 2], 2*Ku1_atom/magnumaf.Constants.mu0/p )
    self.assertAlmostEqual(np_heff[1, 0, 0, 0], 0 )
    self.assertAlmostEqual(np_heff[1, 0, 0, 1], 0 )
    self.assertAlmostEqual(np_heff[1, 0, 0, 2], 0 )

  def test_atomistic_anisotropy_1_2_1_z_z(self):
    mesh=magnumaf.Mesh(1, 2, 1, self.dx, self.dx, self.dx)
    p =self.p
    Ku1_atom =self.k
    m=af.constant(0.0, 1, 2, 1, 3, dtype=af.Dtype.f64)
    m[0, 0, 0, 0] = 0
    m[0, 0, 0, 1] = 0
    m[0, 0, 0, 2] = 1

    m[0, 1, 0, 0] = 0
    m[0, 1, 0, 1] = 0
    m[0, 1, 0, 2] = 1

    state=magnumaf.State(mesh, Ms = p, m = m)
    atom_ani=magnumaf.AtomisticUniaxialAnisotropyField(Ku1_atom, [0, 0, 1])

    self.assertAlmostEqual(atom_ani.Energy_in_J(state), -2*Ku1_atom)

    af_heff = atom_ani.h(state)
    np_heff = af_heff.__array__()

    self.assertAlmostEqual(np_heff[0, 0, 0, 0], 0 )
    self.assertAlmostEqual(np_heff[0, 0, 0, 1], 0 )
    self.assertAlmostEqual(np_heff[0, 0, 0, 2], 2*Ku1_atom/magnumaf.Constants.mu0/p )
    self.assertAlmostEqual(np_heff[0, 1, 0, 0], 0 )
    self.assertAlmostEqual(np_heff[0, 1, 0, 1], 0 )
    self.assertAlmostEqual(np_heff[0, 1, 0, 2], 2*Ku1_atom/magnumaf.Constants.mu0/p )

if __name__ == '__main__':
  unittest.main()
