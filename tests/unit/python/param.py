import unittest
import arrayfire as af
import magnumaf
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
    material=magnumaf.Material(alpha=self.alpha, T=self.T, ms=self.ms, A=self.A, D=self.D, Ku1=self.Ku1, J_atom=self.J_atom, D_atom=self.D_atom, Ku1_axis=self.Ku1_axis, Ku1_atom_axis=self.Ku1_atom_axis, D_axis=self.D_axis, D_atom_axis=self.D_atom_axis)

    self.assertEqual(self.alpha, material.alpha)
    self.assertEqual(self.T, material.T)
    self.assertEqual(self.ms, material.ms)
    self.assertEqual(self.A, material.A)
    self.assertEqual(self.D, material.D)
    self.assertEqual(self.Ku1, material.Ku1)
    self.assertEqual(self.J_atom, material.J_atom)
    self.assertEqual(self.D_atom, material.D_atom)
    self.assertEqual((1.,0.,0.), material.Ku1_axis)
    self.assertEqual((0.,1.,0.), material.Ku1_atom_axis)
    self.assertEqual((1./sqrt(2.),1./sqrt(2.),0.), material.D_axis)
    self.assertEqual((1./sqrt(3.),1./sqrt(3.),1./sqrt(3.)), material.D_atom_axis)

  def test_param_default_initialization(self):
    material=magnumaf.Material()
    self.assertEqual(4e-7*pi, magnumaf.Constants.mu0           )
    self.assertEqual(221276.14886372554, magnumaf.Constants.gamma         )
    self.assertEqual((0.,0.,1.), material.Ku1_axis)
    self.assertEqual(0., material.ms            )
    self.assertEqual(0., material.alpha         )
    self.assertEqual(0., material.A             )
    self.assertEqual(0., material.p             )
    self.assertEqual(0., material.D             )
    self.assertEqual(0., material.Ku1           )
    self.assertEqual(0., material.J_atom        )
    self.assertEqual(0., material.D_atom        )
    self.assertEqual(0., material.Ku1_atom        )

if __name__ == '__main__':
  unittest.main()
