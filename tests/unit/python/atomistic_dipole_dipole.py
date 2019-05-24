import unittest
import arrayfire as af
import numpy as np
import magnumaf
import math

class AtomisticDipoleDipoleTest(unittest.TestCase):
  # arbitrary parameters:
  p  = 9.3e-24
  dx = 2.7e-10

  def test_atomistic_dipole_dipole_2_1_1_z_z(self):
    mesh=magnumaf.Mesh(2, 1, 1, self.dx, self.dx, self.dx)
    material=magnumaf.Material()
    material.p =self.p
    m=af.constant(0.0,2,1,1,3,dtype=af.Dtype.f64)
    m[0,0,0,0] = 0
    m[0,0,0,1] = 0
    m[0,0,0,2] = 1

    m[1,0,0,0] = 0
    m[1,0,0,1] = 0
    m[1,0,0,2] = 1

    pystate=magnumaf.State(mesh, Ms = 0, m = m, material = material)
    atom_demag=magnumaf.AtomisticDipoleDipoleField(mesh)
    Llg=magnumaf.LLGIntegrator(alpha = 1, terms = [atom_demag])

    self.assertEqual(Llg.E(pystate), material.p**2 * magnumaf.Constants.mu0/(4.*math.pi)/self.dx**3)

    af_heff = Llg.h(pystate)
    np_heff = af_heff.__array__()

    self.assertAlmostEqual(np_heff[0,0,0,0], 0 )
    self.assertAlmostEqual(np_heff[0,0,0,1], 0 )
    self.assertAlmostEqual(np_heff[0,0,0,2], -material.p /4. /math.pi /self.dx**3 )
    self.assertAlmostEqual(np_heff[1,0,0,0], 0 )
    self.assertAlmostEqual(np_heff[1,0,0,1], 0 )
    self.assertAlmostEqual(np_heff[1,0,0,2], -material.p /4. /math.pi /self.dx**3 )

  def test_atomistic_dipole_dipole_2_1_1_z_x(self):
    mesh=magnumaf.Mesh(2, 1, 1, self.dx, self.dx, self.dx)
    material=magnumaf.Material()
    material.p =self.p
    m=af.constant(0.0,2,1,1,3,dtype=af.Dtype.f64)
    m[0,0,0,0] = 0
    m[0,0,0,1] = 0
    m[0,0,0,2] = 1

    m[1,0,0,0] = 1
    m[1,0,0,1] = 0
    m[1,0,0,2] = 0

    pystate=magnumaf.State(mesh, Ms = 0, m = m, material = material)
    atom_demag=magnumaf.AtomisticDipoleDipoleField(mesh)
    Llg=magnumaf.LLGIntegrator(alpha = 1, terms = [atom_demag])

    self.assertAlmostEqual(Llg.E(pystate), 0)

    af_heff = Llg.h(pystate)
    np_heff = af_heff.__array__()

    self.assertAlmostEqual(np_heff[0,0,0,0], 2*material.p /4. /math.pi /self.dx**3  )
    self.assertAlmostEqual(np_heff[0,0,0,1], 0 )
    self.assertAlmostEqual(np_heff[0,0,0,2], 0 )
    self.assertAlmostEqual(np_heff[1,0,0,0], 0 )
    self.assertAlmostEqual(np_heff[1,0,0,1], 0 )
    self.assertAlmostEqual(np_heff[1,0,0,2], -material.p /4. /math.pi /self.dx**3 )

  def test_atomistic_dipole_dipole_2_1_1_z_minus_z(self):
    mesh=magnumaf.Mesh(2, 1, 1, self.dx, self.dx, self.dx)
    material=magnumaf.Material()
    material.p = self.p
    m=af.constant(0.0,2,1,1,3,dtype=af.Dtype.f64)
    m[0,0,0,0] = 0
    m[0,0,0,1] = 0
    m[0,0,0,2] = 1

    m[1,0,0,0] = 0
    m[1,0,0,1] = 0
    m[1,0,0,2] = -1

    pystate=magnumaf.State(mesh, Ms = 0, m = m, material = material)
    atom_demag=magnumaf.AtomisticDipoleDipoleField(mesh)
    Llg=magnumaf.LLGIntegrator(alpha = 1, terms = [atom_demag])

    self.assertEqual(Llg.E(pystate), -material.p**2 * magnumaf.Constants.mu0/(4.*math.pi)/self.dx**3)

    af_heff = Llg.h(pystate)
    np_heff = af_heff.__array__()

    self.assertAlmostEqual(np_heff[0,0,0,0], 0 )
    self.assertAlmostEqual(np_heff[0,0,0,1], 0 )
    self.assertAlmostEqual(np_heff[0,0,0,2], material.p /4. /math.pi /self.dx**3 )
    self.assertAlmostEqual(np_heff[1,0,0,0], 0 )
    self.assertAlmostEqual(np_heff[1,0,0,1], 0 )
    self.assertAlmostEqual(np_heff[1,0,0,2], -material.p /4. /math.pi /self.dx**3 )

  def test_atomistic_dipole_dipole_1_2_1_z_z(self):
    mesh=magnumaf.Mesh(1, 2, 1, self.dx, self.dx, self.dx)
    material=magnumaf.Material()
    material.p = self.p
    m=af.constant(0.0,1,2,1,3,dtype=af.Dtype.f64)
    m[0,0,0,0] = 0
    m[0,0,0,1] = 0
    m[0,0,0,2] = 1

    m[0,1,0,0] = 0
    m[0,1,0,1] = 0
    m[0,1,0,2] = 1

    pystate=magnumaf.State(mesh, Ms = 0, m = m, material = material)
    atom_demag=magnumaf.AtomisticDipoleDipoleField(mesh)
    Llg=magnumaf.LLGIntegrator(alpha = 1, terms = [atom_demag])
    
    self.assertEqual(Llg.E(pystate), material.p**2 * magnumaf.Constants.mu0/(4.*math.pi)/self.dx**3)

  def test_atomistic_dipole_dipole_1_2_1_x_z(self):
    mesh=magnumaf.Mesh(1, 2, 1, self.dx, self.dx, self.dx)
    material=magnumaf.Material()
    material.p =self.p
    m=af.constant(0.0,1,2,1,3,dtype=af.Dtype.f64)
    m[0,0,0,0] = 1
    m[0,0,0,1] = 0
    m[0,0,0,2] = 0

    m[0,1,0,0] = 0
    m[0,1,0,1] = 0
    m[0,1,0,2] = 1

    pystate=magnumaf.State(mesh, Ms = 0, m = m, material = material)
    atom_demag=magnumaf.AtomisticDipoleDipoleField(mesh)
    Llg=magnumaf.LLGIntegrator(alpha = 1, terms = [atom_demag])

    self.assertAlmostEqual(Llg.E(pystate), 0)

  def test_atomistic_dipole_dipole_2_1_1_x_x(self):
    mesh=magnumaf.Mesh(2, 1, 1, self.dx, self.dx, self.dx)
    material=magnumaf.Material()
    material.p =self.p
    m=af.constant(0.0,2,1,1,3,dtype=af.Dtype.f64)
    m[0,0,0,0] = 1
    m[0,0,0,1] = 0
    m[0,0,0,2] = 0

    m[1,0,0,0] = 1
    m[1,0,0,1] = 0
    m[1,0,0,2] = 0

    pystate=magnumaf.State(mesh, Ms = 0, m = m, material = material)
    atom_demag=magnumaf.AtomisticDipoleDipoleField(mesh)
    Llg=magnumaf.LLGIntegrator(alpha = 1, terms = [atom_demag])

    self.assertAlmostEqual(Llg.E(pystate), -material.p**2 * magnumaf.Constants.mu0 /(2.*math.pi) /self.dx**3)

    af_heff = Llg.h(pystate)
    np_heff = af_heff.__array__()

    self.assertAlmostEqual(np_heff[0,0,0,0], material.p/2./math.pi /self.dx**3 )
    self.assertAlmostEqual(np_heff[0,0,0,1], 0 )
    self.assertAlmostEqual(np_heff[0,0,0,2], 0 )
    self.assertAlmostEqual(np_heff[1,0,0,0], material.p/2./math.pi /self.dx**3 )
    self.assertAlmostEqual(np_heff[1,0,0,1], 0 )
    self.assertAlmostEqual(np_heff[1,0,0,2], 0 )

if __name__ == '__main__':
  unittest.main()
