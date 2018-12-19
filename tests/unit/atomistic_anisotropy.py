import unittest
import arrayfire as af
import numpy as np
import magnum_af
import math

class AtomisticAnisotropyTest(unittest.TestCase):
  # arbitrary parameters:
  p  = 1e-20
  k  = 1e-20
  dx = 2.715e-10

  def test_atomistic_anisotropy_2_1_1_z_z(self):
    mesh=magnum_af.pyMesh(2, 1, 1, self.dx, self.dx, self.dx)
    param=magnum_af.pyParam()
    param.set_p(self.p)
    param.set_Ku1_atom(self.k)
    m=af.constant(0.0,2,1,1,3,dtype=af.Dtype.f64)
    m[0,0,0,0] = 0
    m[0,0,0,1] = 0
    m[0,0,0,2] = 1

    m[1,0,0,0] = 0
    m[1,0,0,1] = 0
    m[1,0,0,2] = 1

    pystate=magnum_af.pyState(mesh,param,m)
    atom_ani=magnum_af.pyATOMISTIC_ANISOTROPY(mesh, param)
    Llg=magnum_af.pyLLG(atom_ani)

    self.assertAlmostEqual(Llg.get_E(pystate), -2*param.get_Ku1_atom())

    af_heff = Llg.get_fheff(pystate)
    np_heff = af_heff.__array__()

    self.assertAlmostEqual(np_heff[0,0,0,0], 0 )
    self.assertAlmostEqual(np_heff[0,0,0,1], 0 )
    self.assertAlmostEqual(np_heff[0,0,0,2], 2*param.get_Ku1_atom()/param.get_mu0()/param.get_p() )
    self.assertAlmostEqual(np_heff[1,0,0,0], 0 )
    self.assertAlmostEqual(np_heff[1,0,0,1], 0 )
    self.assertAlmostEqual(np_heff[1,0,0,2], 2*param.get_Ku1_atom()/param.get_mu0()/param.get_p() )

  def test_atomistic_anisotropy_2_1_1_z_x(self):
    mesh=magnum_af.pyMesh(2, 1, 1, self.dx, self.dx, self.dx)
    param=magnum_af.pyParam()
    param.set_p(self.p)
    param.set_Ku1_atom(self.k)
    m=af.constant(0.0,2,1,1,3,dtype=af.Dtype.f64)
    m[0,0,0,0] = 0
    m[0,0,0,1] = 0
    m[0,0,0,2] = 1

    m[1,0,0,0] = 1
    m[1,0,0,1] = 0
    m[1,0,0,2] = 0

    pystate=magnum_af.pyState(mesh,param,m)
    atom_ani=magnum_af.pyATOMISTIC_ANISOTROPY(mesh, param)
    Llg=magnum_af.pyLLG(atom_ani)

    self.assertAlmostEqual(Llg.get_E(pystate), -param.get_Ku1_atom())

    af_heff = Llg.get_fheff(pystate)
    np_heff = af_heff.__array__()

    self.assertAlmostEqual(np_heff[0,0,0,0], 0 )
    self.assertAlmostEqual(np_heff[0,0,0,1], 0 )
    self.assertAlmostEqual(np_heff[0,0,0,2], 2*param.get_Ku1_atom()/param.get_mu0()/param.get_p() )
    self.assertAlmostEqual(np_heff[1,0,0,0], 0 )
    self.assertAlmostEqual(np_heff[1,0,0,1], 0 )
    self.assertAlmostEqual(np_heff[1,0,0,2], 0 )

  def test_atomistic_anisotropy_1_2_1_z_z(self):
    mesh=magnum_af.pyMesh(1, 2, 1, self.dx, self.dx, self.dx)
    param=magnum_af.pyParam()
    param.set_p(self.p)
    param.set_Ku1_atom(self.k)
    m=af.constant(0.0,1,2,1,3,dtype=af.Dtype.f64)
    m[0,0,0,0] = 0
    m[0,0,0,1] = 0
    m[0,0,0,2] = 1

    m[0,1,0,0] = 0
    m[0,1,0,1] = 0
    m[0,1,0,2] = 1

    pystate=magnum_af.pyState(mesh,param,m)
    atom_ani=magnum_af.pyATOMISTIC_ANISOTROPY(mesh, param)
    Llg=magnum_af.pyLLG(atom_ani)

    self.assertAlmostEqual(Llg.get_E(pystate), -2*param.get_Ku1_atom())

    af_heff = Llg.get_fheff(pystate)
    np_heff = af_heff.__array__()

    self.assertAlmostEqual(np_heff[0,0,0,0], 0 )
    self.assertAlmostEqual(np_heff[0,0,0,1], 0 )
    self.assertAlmostEqual(np_heff[0,0,0,2], 2*param.get_Ku1_atom()/param.get_mu0()/param.get_p() )
    self.assertAlmostEqual(np_heff[0,1,0,0], 0 )
    self.assertAlmostEqual(np_heff[0,1,0,1], 0 )
    self.assertAlmostEqual(np_heff[0,1,0,2], 2*param.get_Ku1_atom()/param.get_mu0()/param.get_p() )

if __name__ == '__main__':
  unittest.main()
