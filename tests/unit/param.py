import unittest
import arrayfire as af
import magnum_af
from math import pi, sqrt

class ParamTest(unittest.TestCase):
  alpha  = 1.5
  T      = 2.5
  ms     = 3.5
  A      = 4.5
  D      = 5.5
  Ku1    = 6.5
  J_atom = 7.5
  D_atom = 8.5

  Ku1_axis = [2.0, 0., 0.]
  Ku1_atom_axis = [0., 3., 0.]
  D_axis = [1., 1., 0.]
  D_atom_axis = [1., 1., 1.]

  def test_param_initialization(self):
    param=magnum_af.Param(alpha=self.alpha, T=self.T, ms=self.ms, A=self.A, D=self.D, Ku1=self.Ku1, J_atom=self.J_atom, D_atom=self.D_atom, Ku1_axis=self.Ku1_axis, Ku1_atom_axis=self.Ku1_atom_axis, D_axis=self.D_axis, D_atom_axis=self.D_atom_axis)

    self.assertEqual(self.alpha, param.alpha)
    self.assertEqual(self.T, param.T)
    self.assertEqual(self.ms, param.ms)
    self.assertEqual(self.A, param.A)
    self.assertEqual(self.D, param.D)
    self.assertEqual(self.Ku1, param.Ku1)
    self.assertEqual(self.J_atom, param.J_atom)
    self.assertEqual(self.D_atom, param.D_atom)
    self.assertEqual((1.,0.,0.), param.Ku1_axis)
    self.assertEqual((0.,1.,0.), param.Ku1_atom_axis)
    self.assertEqual((1./sqrt(2.),1./sqrt(2.),0.), param.D_axis)
    self.assertEqual((1./sqrt(3.),1./sqrt(3.),1./sqrt(3.)), param.D_atom_axis)

  def test_param_default_initialization(self):
    param=magnum_af.Param()
    self.assertEqual(4e-7*pi, param.mu0           )
    self.assertEqual(221276.1488637255, param.gamma         )
    self.assertEqual((0.,0.,1.), param.Ku1_axis)
    self.assertEqual(0., param.ms            )
    self.assertEqual(0., param.alpha         )
    self.assertEqual(0., param.A             )
    self.assertEqual(0., param.p             )
    self.assertEqual(0., param.D             )
    self.assertEqual(0., param.Ku1           )
    self.assertEqual(0., param.J_atom        )
    self.assertEqual(0., param.D_atom        )
    self.assertEqual(0., param.Ku1_atom        )

if __name__ == '__main__':
  unittest.main()
