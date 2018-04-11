import unittest
import arrayfire as af
import numpy as np
import pth_mag
import math

class AtomisticAnisotropyTest(unittest.TestCase):
  # arbitrary parameters:
  p  = 1e-20
  k  = 1e-20
  dx = 2.715e-10

  def test_atomistic_dipole_dipole_2_1_1_z_z(self):
    mesh=pth_mag.pyMesh(2, 1, 1, self.dx, self.dx, self.dx)
    param=pth_mag.pyParam()
    param.p (self.p)
    param.K_atom (self.k)
    m=af.constant(0.0,2,1,1,3,dtype=af.Dtype.f64)
    m[0,0,0,0] = 0
    m[0,0,0,1] = 0
    m[0,0,0,2] = 1

    m[1,0,0,0] = 0
    m[1,0,0,1] = 0
    m[1,0,0,2] = 1

    pystate=pth_mag.pyState(mesh,param,m)
    atom_ani=pth_mag.pyATOMISTIC_ANISOTROPY(mesh, param)
    Llg=pth_mag.pyLLG(pystate,atom_ani)

    self.assertAlmostEqual(Llg.print_E(pystate), -2*param.print_K_atom())

    af_heff = Llg.get_fheff(pystate)
    np_heff = af_heff.__array__()

    self.assertAlmostEqual(np_heff[0,0,0,0], 0 )
    self.assertAlmostEqual(np_heff[0,0,0,1], 0 )
    self.assertAlmostEqual(np_heff[0,0,0,2], 2*param.print_K_atom()/param.print_mu0()/param.print_p() )
    self.assertAlmostEqual(np_heff[1,0,0,0], 0 )
    self.assertAlmostEqual(np_heff[1,0,0,1], 0 )
    self.assertAlmostEqual(np_heff[1,0,0,2], 2*param.print_K_atom()/param.print_mu0()/param.print_p() )

  def test_atomistic_dipole_dipole_2_1_1_z_x(self):
    mesh=pth_mag.pyMesh(2, 1, 1, self.dx, self.dx, self.dx)
    param=pth_mag.pyParam()
    param.p (self.p)
    param.K_atom (self.k)
    m=af.constant(0.0,2,1,1,3,dtype=af.Dtype.f64)
    m[0,0,0,0] = 0
    m[0,0,0,1] = 0
    m[0,0,0,2] = 1

    m[1,0,0,0] = 1
    m[1,0,0,1] = 0
    m[1,0,0,2] = 0

    pystate=pth_mag.pyState(mesh,param,m)
    atom_ani=pth_mag.pyATOMISTIC_ANISOTROPY(mesh, param)
    Llg=pth_mag.pyLLG(pystate,atom_ani)

    self.assertAlmostEqual(Llg.print_E(pystate), -param.print_K_atom())

    af_heff = Llg.get_fheff(pystate)
    np_heff = af_heff.__array__()

    self.assertAlmostEqual(np_heff[0,0,0,0], 0 )
    self.assertAlmostEqual(np_heff[0,0,0,1], 0 )
    self.assertAlmostEqual(np_heff[0,0,0,2], 2*param.print_K_atom()/param.print_mu0()/param.print_p() )
    self.assertAlmostEqual(np_heff[1,0,0,0], 0 )
    self.assertAlmostEqual(np_heff[1,0,0,1], 0 )
    self.assertAlmostEqual(np_heff[1,0,0,2], 0 )

  def test_atomistic_dipole_dipole_1_2_1_z_z(self):
    mesh=pth_mag.pyMesh(1, 2, 1, self.dx, self.dx, self.dx)
    param=pth_mag.pyParam()
    param.p (self.p)
    param.K_atom (self.k)
    m=af.constant(0.0,1,2,1,3,dtype=af.Dtype.f64)
    m[0,0,0,0] = 0
    m[0,0,0,1] = 0
    m[0,0,0,2] = 1

    m[0,1,0,0] = 0
    m[0,1,0,1] = 0
    m[0,1,0,2] = 1

    pystate=pth_mag.pyState(mesh,param,m)
    atom_ani=pth_mag.pyATOMISTIC_ANISOTROPY(mesh, param)
    Llg=pth_mag.pyLLG(pystate,atom_ani)

    self.assertAlmostEqual(Llg.print_E(pystate), -2*param.print_K_atom())

    af_heff = Llg.get_fheff(pystate)
    np_heff = af_heff.__array__()

    self.assertAlmostEqual(np_heff[0,0,0,0], 0 )
    self.assertAlmostEqual(np_heff[0,0,0,1], 0 )
    self.assertAlmostEqual(np_heff[0,0,0,2], 2*param.print_K_atom()/param.print_mu0()/param.print_p() )
    self.assertAlmostEqual(np_heff[0,1,0,0], 0 )
    self.assertAlmostEqual(np_heff[0,1,0,1], 0 )
    self.assertAlmostEqual(np_heff[0,1,0,2], 2*param.print_K_atom()/param.print_mu0()/param.print_p() )

if __name__ == '__main__':
  unittest.main()
