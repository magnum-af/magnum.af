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
    p =self.p
    m=af.constant(0.0, 2, 1, 1, 3, dtype=af.Dtype.f64)
    m[0, 0, 0, 0] = 0
    m[0, 0, 0, 1] = 0
    m[0, 0, 0, 2] = 1

    m[1, 0, 0, 0] = 0
    m[1, 0, 0, 1] = 0
    m[1, 0, 0, 2] = 1

    state=magnumaf.State(mesh, Ms = p, m = m)
    atom_demag=magnumaf.AtomisticDipoleDipoleField(mesh)

    self.assertEqual(atom_demag.Energy_in_J(state), p**2 * magnumaf.Constants.mu0/(4.*math.pi)/self.dx**3)

    af_heff = atom_demag.H_in_Apm(state)
    np_heff = af_heff.__array__()

    self.assertAlmostEqual(np_heff[0, 0, 0, 0], 0 )
    self.assertAlmostEqual(np_heff[0, 0, 0, 1], 0 )
    self.assertAlmostEqual(np_heff[0, 0, 0, 2], -p /4. /math.pi /self.dx**3 )
    self.assertAlmostEqual(np_heff[1, 0, 0, 0], 0 )
    self.assertAlmostEqual(np_heff[1, 0, 0, 1], 0 )
    self.assertAlmostEqual(np_heff[1, 0, 0, 2], -p /4. /math.pi /self.dx**3 )

  def test_atomistic_dipole_dipole_2_1_1_z_x(self):
    mesh=magnumaf.Mesh(2, 1, 1, self.dx, self.dx, self.dx)
    p =self.p
    m=af.constant(0.0, 2, 1, 1, 3, dtype=af.Dtype.f64)
    m[0, 0, 0, 0] = 0
    m[0, 0, 0, 1] = 0
    m[0, 0, 0, 2] = 1

    m[1, 0, 0, 0] = 1
    m[1, 0, 0, 1] = 0
    m[1, 0, 0, 2] = 0

    state=magnumaf.State(mesh, Ms = p, m = m)
    atom_demag=magnumaf.AtomisticDipoleDipoleField(mesh)

    self.assertAlmostEqual(atom_demag.Energy_in_J(state), 0)

    af_heff = atom_demag.H_in_Apm(state)
    np_heff = af_heff.__array__()

    self.assertAlmostEqual(np_heff[0, 0, 0, 0], 2*p /4. /math.pi /self.dx**3  )
    self.assertAlmostEqual(np_heff[0, 0, 0, 1], 0 )
    self.assertAlmostEqual(np_heff[0, 0, 0, 2], 0 )
    self.assertAlmostEqual(np_heff[1, 0, 0, 0], 0 )
    self.assertAlmostEqual(np_heff[1, 0, 0, 1], 0 )
    self.assertAlmostEqual(np_heff[1, 0, 0, 2], -p /4. /math.pi /self.dx**3 )

  def test_atomistic_dipole_dipole_2_1_1_z_minus_z(self):
    mesh=magnumaf.Mesh(2, 1, 1, self.dx, self.dx, self.dx)
    p = self.p
    m=af.constant(0.0, 2, 1, 1, 3, dtype=af.Dtype.f64)
    m[0, 0, 0, 0] = 0
    m[0, 0, 0, 1] = 0
    m[0, 0, 0, 2] = 1

    m[1, 0, 0, 0] = 0
    m[1, 0, 0, 1] = 0
    m[1, 0, 0, 2] = -1

    state=magnumaf.State(mesh, Ms = p, m = m)
    atom_demag=magnumaf.AtomisticDipoleDipoleField(mesh)

    self.assertEqual(atom_demag.Energy_in_J(state), -p**2 * magnumaf.Constants.mu0/(4.*math.pi)/self.dx**3)

    af_heff = atom_demag.H_in_Apm(state)
    np_heff = af_heff.__array__()

    self.assertAlmostEqual(np_heff[0, 0, 0, 0], 0 )
    self.assertAlmostEqual(np_heff[0, 0, 0, 1], 0 )
    self.assertAlmostEqual(np_heff[0, 0, 0, 2], p /4. /math.pi /self.dx**3 )
    self.assertAlmostEqual(np_heff[1, 0, 0, 0], 0 )
    self.assertAlmostEqual(np_heff[1, 0, 0, 1], 0 )
    self.assertAlmostEqual(np_heff[1, 0, 0, 2], -p /4. /math.pi /self.dx**3 )

  def test_atomistic_dipole_dipole_1_2_1_z_z(self):
    mesh=magnumaf.Mesh(1, 2, 1, self.dx, self.dx, self.dx)
    p = self.p
    m=af.constant(0.0, 1, 2, 1, 3, dtype=af.Dtype.f64)
    m[0, 0, 0, 0] = 0
    m[0, 0, 0, 1] = 0
    m[0, 0, 0, 2] = 1

    m[0, 1, 0, 0] = 0
    m[0, 1, 0, 1] = 0
    m[0, 1, 0, 2] = 1

    state=magnumaf.State(mesh, Ms = p, m = m)
    atom_demag=magnumaf.AtomisticDipoleDipoleField(mesh)

    self.assertEqual(atom_demag.Energy_in_J(state), p**2 * magnumaf.Constants.mu0/(4.*math.pi)/self.dx**3)

  def test_atomistic_dipole_dipole_1_2_1_x_z(self):
    mesh=magnumaf.Mesh(1, 2, 1, self.dx, self.dx, self.dx)
    p =self.p
    m=af.constant(0.0, 1, 2, 1, 3, dtype=af.Dtype.f64)
    m[0, 0, 0, 0] = 1
    m[0, 0, 0, 1] = 0
    m[0, 0, 0, 2] = 0

    m[0, 1, 0, 0] = 0
    m[0, 1, 0, 1] = 0
    m[0, 1, 0, 2] = 1

    state=magnumaf.State(mesh, Ms = p, m = m)
    atom_demag=magnumaf.AtomisticDipoleDipoleField(mesh)

    self.assertAlmostEqual(atom_demag.Energy_in_J(state), 0)

  def test_atomistic_dipole_dipole_2_1_1_x_x(self):
    mesh=magnumaf.Mesh(2, 1, 1, self.dx, self.dx, self.dx)
    p =self.p
    m=af.constant(0.0, 2, 1, 1, 3, dtype=af.Dtype.f64)
    m[0, 0, 0, 0] = 1
    m[0, 0, 0, 1] = 0
    m[0, 0, 0, 2] = 0

    m[1, 0, 0, 0] = 1
    m[1, 0, 0, 1] = 0
    m[1, 0, 0, 2] = 0

    state=magnumaf.State(mesh, Ms = p, m = m)
    atom_demag=magnumaf.AtomisticDipoleDipoleField(mesh)

    self.assertAlmostEqual(atom_demag.Energy_in_J(state), -p**2 * magnumaf.Constants.mu0 /(2.*math.pi) /self.dx**3)

    af_heff = atom_demag.H_in_Apm(state)
    np_heff = af_heff.__array__()

    self.assertAlmostEqual(np_heff[0, 0, 0, 0], p/2./math.pi /self.dx**3 )
    self.assertAlmostEqual(np_heff[0, 0, 0, 1], 0 )
    self.assertAlmostEqual(np_heff[0, 0, 0, 2], 0 )
    self.assertAlmostEqual(np_heff[1, 0, 0, 0], p/2./math.pi /self.dx**3 )
    self.assertAlmostEqual(np_heff[1, 0, 0, 1], 0 )
    self.assertAlmostEqual(np_heff[1, 0, 0, 2], 0 )

if __name__ == '__main__':
  unittest.main()
