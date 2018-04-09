import unittest
import arrayfire as af
import pth_mag
import math

class AtomisticAnisotropyTest(unittest.TestCase):
    #TODO also test H field
  p=9.274009994e-24
  dx = 2.715e-10

  def test_atomistic_anisotropy_2_1_1_z_z(self):
    mesh=pth_mag.pyMesh(2, 1, 1, self.dx, self.dx, self.dx)
    param=pth_mag.pyParam()
    param.alpha (1)
    param.p (self.p)
    m=af.constant(0.0,2,1,1,3,dtype=af.Dtype.f64)
    m[0,0,0,0] = 0
    m[0,0,0,1] = 0
    m[0,0,0,2] = 1

    m[1,0,0,0] = 0
    m[1,0,0,1] = 0
    m[1,0,0,2] = 1

    pystate=pth_mag.pyState(mesh,param,m)
    atom_demag=pth_mag.pyATOMISTIC_DEMAG(mesh)
    Llg=pth_mag.pyLLG(pystate,atom_demag)
    analytical = -param.print_p()**2 * param.print_mu0()/(4.*math.pi)/self.dx**3

    self.assertEqual(Llg.print_E(pystate),analytical)

  def test_atomistic_anisotropy_2_1_1_z_x(self):
    mesh=pth_mag.pyMesh(2, 1, 1, self.dx, self.dx, self.dx)
    param=pth_mag.pyParam()
    param.alpha (1)
    param.p (self.p)
    m=af.constant(0.0,2,1,1,3,dtype=af.Dtype.f64)
    m[0,0,0,0] = 0
    m[0,0,0,1] = 0
    m[0,0,0,2] = 1

    m[1,0,0,0] = 1
    m[1,0,0,1] = 0
    m[1,0,0,2] = 0

    pystate=pth_mag.pyState(mesh,param,m)
    atom_demag=pth_mag.pyATOMISTIC_DEMAG(mesh)
    Llg=pth_mag.pyLLG(pystate,atom_demag)
    analytical = 0

    self.assertAlmostEqual(Llg.print_E(pystate),analytical)

  def test_atomistic_anisotropy_2_1_1_z_minus_z(self):
    mesh=pth_mag.pyMesh(2, 1, 1, self.dx, self.dx, self.dx)
    param=pth_mag.pyParam()
    param.alpha (1)
    param.p (self.p)
    m=af.constant(0.0,2,1,1,3,dtype=af.Dtype.f64)
    m[0,0,0,0] = 0
    m[0,0,0,1] = 0
    m[0,0,0,2] = 1

    m[1,0,0,0] = 0
    m[1,0,0,1] = 0
    m[1,0,0,2] = -1

    pystate=pth_mag.pyState(mesh,param,m)
    atom_demag=pth_mag.pyATOMISTIC_DEMAG(mesh)
    Llg=pth_mag.pyLLG(pystate,atom_demag)
    analytical = param.print_p()**2 * param.print_mu0()/(4.*math.pi)/self.dx**3

    self.assertEqual(Llg.print_E(pystate),analytical)

  def test_atomistic_anisotropy_1_2_1_z_z(self):
    mesh=pth_mag.pyMesh(2, 1, 1, self.dx, self.dx, self.dx)
    param=pth_mag.pyParam()
    param.alpha (1)
    param.p (self.p)
    m=af.constant(0.0,2,1,1,3,dtype=af.Dtype.f64)
    m[0,0,0,0] = 0
    m[0,0,0,1] = 0
    m[0,0,0,2] = 1

    m[1,0,0,0] = 0
    m[1,0,0,1] = 0
    m[1,0,0,2] = 1

    pystate=pth_mag.pyState(mesh,param,m)
    atom_demag=pth_mag.pyATOMISTIC_DEMAG(mesh)
    Llg=pth_mag.pyLLG(pystate,atom_demag)
    analytical = -param.print_p()**2 * param.print_mu0()/(4.*math.pi)/self.dx**3

    self.assertEqual(Llg.print_E(pystate),analytical)

  def test_atomistic_anisotropy_1_2_1_x_z(self):
    mesh=pth_mag.pyMesh(2, 1, 1, self.dx, self.dx, self.dx)
    param=pth_mag.pyParam()
    param.alpha (1)
    param.p (self.p)
    m=af.constant(0.0,2,1,1,3,dtype=af.Dtype.f64)
    m[0,0,0,0] = 1
    m[0,0,0,1] = 0
    m[0,0,0,2] = 0

    m[1,0,0,0] = 0
    m[1,0,0,1] = 0
    m[1,0,0,2] = 1

    pystate=pth_mag.pyState(mesh,param,m)
    atom_demag=pth_mag.pyATOMISTIC_DEMAG(mesh)
    Llg=pth_mag.pyLLG(pystate,atom_demag)
    analytical = 0

    self.assertAlmostEqual(Llg.print_E(pystate),analytical)

if __name__ == '__main__':
  unittest.main()
