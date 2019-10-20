import unittest
import arrayfire as af
import numpy as np
import magnumaf
import math

class AtomisticDmiFieldTest(unittest.TestCase):
  # arbitrary parameters:
  p  = 1e-20
  D  = 1e20
  dx = 1.1

  def test_atomistic_dmi_2_1_1_z_z(self):
    mesh=magnumaf.Mesh(2, 1, 1, self.dx, self.dx, self.dx)
    material=magnumaf.Material()
    material.p =self.p
    material.D_atom =self.D
    m=af.constant(0.0, 2, 1, 1, 3, dtype=af.Dtype.f32)
    m[0, 0, 0, 0] = 0
    m[0, 0, 0, 1] = 0
    m[0, 0, 0, 2] = 1

    m[1, 0, 0, 0] = 0
    m[1, 0, 0, 1] = 0
    m[1, 0, 0, 2] = 1

    state=magnumaf.State(mesh, Ms = 0, m = m, material = material)
    atom_dmi=magnumaf.AtomisticDmiField(mesh, material)
    Llg=magnumaf.LLGIntegrator(alpha = 0, terms = [atom_dmi])

    self.assertAlmostEqual(Llg.E(state), 0)

    af_heff = Llg.h(state)
    np_heff = af_heff.__array__()

    self.assertAlmostEqual(np_heff[0, 0, 0, 0], -material.D_atom/magnumaf.Constants.mu0/material.p )
    self.assertLess(math.fabs((np_heff[0, 0, 0, 0] - (-material.D_atom/magnumaf.Constants.mu0/material.p))/np_heff[0, 0, 0, 0]), 1e-15 )
    self.assertAlmostEqual(np_heff[0, 0, 0, 1], 0 )
    self.assertAlmostEqual(np_heff[0, 0, 0, 2], 0 )
    self.assertAlmostEqual(np_heff[1, 0, 0, 0], material.D_atom/magnumaf.Constants.mu0/material.p)
    self.assertLess(math.fabs((np_heff[1, 0, 0, 0] - material.D_atom/magnumaf.Constants.mu0/material.p)/np_heff[1, 0, 0, 0]), 1e-15)
    self.assertAlmostEqual(np_heff[1, 0, 0, 1], 0 )
    self.assertAlmostEqual(np_heff[1, 0, 0, 2], 0 )

  def test_atomistic_dmi_2_1_1_z_x(self):
    mesh=magnumaf.Mesh(2, 1, 1, self.dx, self.dx, self.dx)
    material=magnumaf.Material()
    material.p =self.p
    material.D_atom =self.D
    m=af.constant(0.0, 2, 1, 1, 3, dtype=af.Dtype.f32)
    m[0, 0, 0, 0] = 0
    m[0, 0, 0, 1] = 0
    m[0, 0, 0, 2] = 1

    m[1, 0, 0, 0] = 1
    m[1, 0, 0, 1] = 0
    m[1, 0, 0, 2] = 0

    state=magnumaf.State(mesh, Ms = 0, m = m, material = material)
    atom_dmi=magnumaf.AtomisticDmiField(mesh, material)
    Llg=magnumaf.LLGIntegrator(alpha = 0, terms = [atom_dmi])

    self.assertAlmostEqual(Llg.E(state), - material.D_atom)

    af_heff = Llg.h(state)
    np_heff = af_heff.__array__()

    self.assertAlmostEqual(np_heff[0, 0, 0, 0], 0 )
    self.assertAlmostEqual(np_heff[0, 0, 0, 1], 0 )
    self.assertAlmostEqual(np_heff[0, 0, 0, 2], material.D_atom/magnumaf.Constants.mu0/material.p )
    self.assertLess(math.fabs((np_heff[0, 0, 0, 2] - material.D_atom/magnumaf.Constants.mu0/material.p)/np_heff[0, 0, 0, 2]), 1e-15 )

    self.assertAlmostEqual(np_heff[1, 0, 0, 0], material.D_atom/magnumaf.Constants.mu0/material.p)
    self.assertLess(math.fabs((np_heff[1, 0, 0, 0] - material.D_atom/magnumaf.Constants.mu0/material.p)/np_heff[1, 0, 0, 0]), 1e-15)
    self.assertAlmostEqual(np_heff[1, 0, 0, 1], 0 )
    self.assertAlmostEqual(np_heff[1, 0, 0, 2], 0 )

  def test_atomistic_dmi_1_1_2_z_z(self):
    mesh=magnumaf.Mesh(1, 1, 2, self.dx, self.dx, self.dx)
    material=magnumaf.Material()
    material.p =self.p
    material.D_atom =self.D
    m=af.constant(0.0, 1, 1, 2, 3, dtype=af.Dtype.f32)
    m[0, 0, 0, 0] = 0
    m[0, 0, 0, 1] = 0
    m[0, 0, 0, 2] = 1

    m[0, 0, 1, 0] = 1
    m[0, 0, 1, 1] = 0
    m[0, 0, 1, 2] = 0

    state=magnumaf.State(mesh, Ms = 0, m = m, material = material)
    atom_dmi=magnumaf.AtomisticDmiField(mesh, material)
    Llg=magnumaf.LLGIntegrator(alpha = 0, terms = [atom_dmi])

    self.assertAlmostEqual(Llg.E(state), 0)

    af_heff = Llg.h(state)
    np_heff = af_heff.__array__()

    self.assertAlmostEqual(np_heff[0, 0, 0, 0], 0 )
    self.assertAlmostEqual(np_heff[0, 0, 0, 1], 0 )
    self.assertAlmostEqual(np_heff[0, 0, 0, 2], 0 )
    self.assertAlmostEqual(np_heff[0, 0, 1, 0], 0 )
    self.assertAlmostEqual(np_heff[0, 0, 1, 1], 0 )
    self.assertAlmostEqual(np_heff[0, 0, 1, 2], 0 )

if __name__ == '__main__':
  unittest.main()
