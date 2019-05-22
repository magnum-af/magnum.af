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

  def test_atomistic_anisotropy_2_1_1_z_z(self):
    mesh=magnumaf.Mesh(2, 1, 1, self.dx, self.dx, self.dx)
    material=magnumaf.Material()
    material.p=self.p
    material.Ku1_atom=self.k
    m=af.constant(0.0,2,1,1,3,dtype=af.Dtype.f64)
    m[0,0,0,0] = 0
    m[0,0,0,1] = 0
    m[0,0,0,2] = 1

    m[1,0,0,0] = 0
    m[1,0,0,1] = 0
    m[1,0,0,2] = 1

    pystate=magnumaf.State(mesh, Ms = 0, m = m, material = material)
    atom_ani=magnumaf.AtomisticUniaxialAnisotropyField(mesh, material)
    Llg=magnumaf.LLGIntegrator(alpha = 0, terms = [atom_ani])

    self.assertAlmostEqual(Llg.get_E(pystate), -2*material.Ku1_atom)

    af_heff = Llg.get_fheff(pystate)
    np_heff = af_heff.__array__()

    self.assertAlmostEqual(np_heff[0,0,0,0], 0 )
    self.assertAlmostEqual(np_heff[0,0,0,1], 0 )
    self.assertAlmostEqual(np_heff[0,0,0,2], 2*material.Ku1_atom/magnumaf.Constants.mu0/material.p )
    self.assertAlmostEqual(np_heff[1,0,0,0], 0 )
    self.assertAlmostEqual(np_heff[1,0,0,1], 0 )
    self.assertAlmostEqual(np_heff[1,0,0,2], 2*material.Ku1_atom/magnumaf.Constants.mu0/material.p )

  def test_atomistic_anisotropy_2_1_1_z_x(self):
    mesh=magnumaf.Mesh(2, 1, 1, self.dx, self.dx, self.dx)
    material=magnumaf.Material()
    material.p =self.p
    material.Ku1_atom =self.k
    m=af.constant(0.0,2,1,1,3,dtype=af.Dtype.f64)
    m[0,0,0,0] = 0
    m[0,0,0,1] = 0
    m[0,0,0,2] = 1

    m[1,0,0,0] = 1
    m[1,0,0,1] = 0
    m[1,0,0,2] = 0

    pystate=magnumaf.State(mesh, Ms = 0, m = m, material = material)
    atom_ani=magnumaf.AtomisticUniaxialAnisotropyField(mesh, material)
    Llg=magnumaf.LLGIntegrator(alpha = 0, terms = [atom_ani])

    self.assertAlmostEqual(Llg.get_E(pystate), -material.Ku1_atom)

    af_heff = Llg.get_fheff(pystate)
    np_heff = af_heff.__array__()

    self.assertAlmostEqual(np_heff[0,0,0,0], 0 )
    self.assertAlmostEqual(np_heff[0,0,0,1], 0 )
    self.assertAlmostEqual(np_heff[0,0,0,2], 2*material.Ku1_atom/magnumaf.Constants.mu0/material.p )
    self.assertAlmostEqual(np_heff[1,0,0,0], 0 )
    self.assertAlmostEqual(np_heff[1,0,0,1], 0 )
    self.assertAlmostEqual(np_heff[1,0,0,2], 0 )

  def test_atomistic_anisotropy_1_2_1_z_z(self):
    mesh=magnumaf.Mesh(1, 2, 1, self.dx, self.dx, self.dx)
    material=magnumaf.Material()
    material.p =self.p
    material.Ku1_atom =self.k
    m=af.constant(0.0,1,2,1,3,dtype=af.Dtype.f64)
    m[0,0,0,0] = 0
    m[0,0,0,1] = 0
    m[0,0,0,2] = 1

    m[0,1,0,0] = 0
    m[0,1,0,1] = 0
    m[0,1,0,2] = 1

    pystate=magnumaf.State(mesh, Ms = 0, m = m, material = material)
    atom_ani=magnumaf.AtomisticUniaxialAnisotropyField(mesh, material)
    Llg=magnumaf.LLGIntegrator(alpha = 0, terms = [atom_ani])

    self.assertAlmostEqual(Llg.get_E(pystate), -2*material.Ku1_atom)

    af_heff = Llg.get_fheff(pystate)
    np_heff = af_heff.__array__()

    self.assertAlmostEqual(np_heff[0,0,0,0], 0 )
    self.assertAlmostEqual(np_heff[0,0,0,1], 0 )
    self.assertAlmostEqual(np_heff[0,0,0,2], 2*material.Ku1_atom/magnumaf.Constants.mu0/material.p )
    self.assertAlmostEqual(np_heff[0,1,0,0], 0 )
    self.assertAlmostEqual(np_heff[0,1,0,1], 0 )
    self.assertAlmostEqual(np_heff[0,1,0,2], 2*material.Ku1_atom/magnumaf.Constants.mu0/material.p )

if __name__ == '__main__':
  unittest.main()
