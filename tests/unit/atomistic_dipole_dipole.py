import unittest
import arrayfire as af
import numpy as np
import magnum_af
import math

class AtomisticDipoleDipoleTest(unittest.TestCase):
  # arbitrary parameters:
  p  = 9.3e-24
  dx = 2.7e-10

  def test_atomistic_dipole_dipole_2_1_1_z_z(self):
    mesh=magnum_af.Mesh(2, 1, 1, self.dx, self.dx, self.dx)
    param=magnum_af.Param()
    param.alpha =1
    param.p =self.p
    m=af.constant(0.0,2,1,1,3,dtype=af.Dtype.f64)
    m[0,0,0,0] = 0
    m[0,0,0,1] = 0
    m[0,0,0,2] = 1

    m[1,0,0,0] = 0
    m[1,0,0,1] = 0
    m[1,0,0,2] = 1

    pystate=magnum_af.State(mesh,param,m)
    atom_demag=magnum_af.pyATOMISTIC_DEMAG(mesh)
    Llg=magnum_af.pyLLG(atom_demag)

    self.assertEqual(Llg.get_E(pystate), param.p**2 * param.mu0/(4.*math.pi)/self.dx**3)

    af_heff = Llg.get_fheff(pystate)
    np_heff = af_heff.__array__()

    self.assertAlmostEqual(np_heff[0,0,0,0], 0 )
    self.assertAlmostEqual(np_heff[0,0,0,1], 0 )
    self.assertAlmostEqual(np_heff[0,0,0,2], -param.p /4. /math.pi /self.dx**3 )
    self.assertAlmostEqual(np_heff[1,0,0,0], 0 )
    self.assertAlmostEqual(np_heff[1,0,0,1], 0 )
    self.assertAlmostEqual(np_heff[1,0,0,2], -param.p /4. /math.pi /self.dx**3 )

  def test_atomistic_dipole_dipole_2_1_1_z_x(self):
    mesh=magnum_af.Mesh(2, 1, 1, self.dx, self.dx, self.dx)
    param=magnum_af.Param()
    param.alpha =1
    param.p =self.p
    m=af.constant(0.0,2,1,1,3,dtype=af.Dtype.f64)
    m[0,0,0,0] = 0
    m[0,0,0,1] = 0
    m[0,0,0,2] = 1

    m[1,0,0,0] = 1
    m[1,0,0,1] = 0
    m[1,0,0,2] = 0

    pystate=magnum_af.State(mesh,param,m)
    atom_demag=magnum_af.pyATOMISTIC_DEMAG(mesh)
    Llg=magnum_af.pyLLG(atom_demag)

    self.assertAlmostEqual(Llg.get_E(pystate), 0)

    af_heff = Llg.get_fheff(pystate)
    np_heff = af_heff.__array__()

    self.assertAlmostEqual(np_heff[0,0,0,0], 2*param.p /4. /math.pi /self.dx**3  )
    self.assertAlmostEqual(np_heff[0,0,0,1], 0 )
    self.assertAlmostEqual(np_heff[0,0,0,2], 0 )
    self.assertAlmostEqual(np_heff[1,0,0,0], 0 )
    self.assertAlmostEqual(np_heff[1,0,0,1], 0 )
    self.assertAlmostEqual(np_heff[1,0,0,2], -param.p /4. /math.pi /self.dx**3 )

  def test_atomistic_dipole_dipole_2_1_1_z_minus_z(self):
    mesh=magnum_af.Mesh(2, 1, 1, self.dx, self.dx, self.dx)
    param=magnum_af.Param()
    param.alpha = 1
    param.p = self.p
    m=af.constant(0.0,2,1,1,3,dtype=af.Dtype.f64)
    m[0,0,0,0] = 0
    m[0,0,0,1] = 0
    m[0,0,0,2] = 1

    m[1,0,0,0] = 0
    m[1,0,0,1] = 0
    m[1,0,0,2] = -1

    pystate=magnum_af.State(mesh,param,m)
    atom_demag=magnum_af.pyATOMISTIC_DEMAG(mesh)
    Llg=magnum_af.pyLLG(atom_demag)

    self.assertEqual(Llg.get_E(pystate), -param.p**2 * param.mu0/(4.*math.pi)/self.dx**3)

    af_heff = Llg.get_fheff(pystate)
    np_heff = af_heff.__array__()

    self.assertAlmostEqual(np_heff[0,0,0,0], 0 )
    self.assertAlmostEqual(np_heff[0,0,0,1], 0 )
    self.assertAlmostEqual(np_heff[0,0,0,2], param.p /4. /math.pi /self.dx**3 )
    self.assertAlmostEqual(np_heff[1,0,0,0], 0 )
    self.assertAlmostEqual(np_heff[1,0,0,1], 0 )
    self.assertAlmostEqual(np_heff[1,0,0,2], -param.p /4. /math.pi /self.dx**3 )

  def test_atomistic_dipole_dipole_1_2_1_z_z(self):
    mesh=magnum_af.Mesh(1, 2, 1, self.dx, self.dx, self.dx)
    param=magnum_af.Param()
    param.alpha = 1
    param.p = self.p
    m=af.constant(0.0,1,2,1,3,dtype=af.Dtype.f64)
    m[0,0,0,0] = 0
    m[0,0,0,1] = 0
    m[0,0,0,2] = 1

    m[0,1,0,0] = 0
    m[0,1,0,1] = 0
    m[0,1,0,2] = 1

    pystate=magnum_af.State(mesh,param,m)
    atom_demag=magnum_af.pyATOMISTIC_DEMAG(mesh)
    Llg=magnum_af.pyLLG(atom_demag)
    
    self.assertEqual(Llg.get_E(pystate), param.p**2 * param.mu0/(4.*math.pi)/self.dx**3)

  def test_atomistic_dipole_dipole_1_2_1_x_z(self):
    mesh=magnum_af.Mesh(1, 2, 1, self.dx, self.dx, self.dx)
    param=magnum_af.Param()
    param.alpha =1
    param.p =self.p
    m=af.constant(0.0,1,2,1,3,dtype=af.Dtype.f64)
    m[0,0,0,0] = 1
    m[0,0,0,1] = 0
    m[0,0,0,2] = 0

    m[0,1,0,0] = 0
    m[0,1,0,1] = 0
    m[0,1,0,2] = 1

    pystate=magnum_af.State(mesh,param,m)
    atom_demag=magnum_af.pyATOMISTIC_DEMAG(mesh)
    Llg=magnum_af.pyLLG(atom_demag)

    self.assertAlmostEqual(Llg.get_E(pystate), 0)

  def test_atomistic_dipole_dipole_2_1_1_x_x(self):
    mesh=magnum_af.Mesh(2, 1, 1, self.dx, self.dx, self.dx)
    param=magnum_af.Param()
    param.alpha =1
    param.p =self.p
    m=af.constant(0.0,2,1,1,3,dtype=af.Dtype.f64)
    m[0,0,0,0] = 1
    m[0,0,0,1] = 0
    m[0,0,0,2] = 0

    m[1,0,0,0] = 1
    m[1,0,0,1] = 0
    m[1,0,0,2] = 0

    pystate=magnum_af.State(mesh,param,m)
    atom_demag=magnum_af.pyATOMISTIC_DEMAG(mesh)
    Llg=magnum_af.pyLLG(atom_demag)

    self.assertAlmostEqual(Llg.get_E(pystate), -param.p**2 * param.mu0 /(2.*math.pi) /self.dx**3)

    af_heff = Llg.get_fheff(pystate)
    np_heff = af_heff.__array__()

    self.assertAlmostEqual(np_heff[0,0,0,0], param.p/2./math.pi /self.dx**3 )
    self.assertAlmostEqual(np_heff[0,0,0,1], 0 )
    self.assertAlmostEqual(np_heff[0,0,0,2], 0 )
    self.assertAlmostEqual(np_heff[1,0,0,0], param.p/2./math.pi /self.dx**3 )
    self.assertAlmostEqual(np_heff[1,0,0,1], 0 )
    self.assertAlmostEqual(np_heff[1,0,0,2], 0 )

if __name__ == '__main__':
  unittest.main()
