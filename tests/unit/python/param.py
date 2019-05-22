import unittest
import arrayfire as af
import magnumaf
from math import pi, sqrt

class ParamTest(unittest.TestCase):
  alpha  = 1.5
  T      = 2.5
  ms     = 3.5
  D      = 5.5
  J_atom = 7.5
  D_atom = 8.5

  Ku1_atom_axis = [0., 3., 0.]
  D_axis = [1., 1., 0.]
  D_atom_axis = [1., 1., 1.]

  def test_param_initialization(self):
    material=magnumaf.Material(D=self.D, J_atom=self.J_atom, D_atom=self.D_atom, Ku1_atom_axis=self.Ku1_atom_axis, D_axis=self.D_axis, D_atom_axis=self.D_atom_axis)

    self.assertEqual(self.D, material.D)
    self.assertEqual(self.J_atom, material.J_atom)
    self.assertEqual(self.D_atom, material.D_atom)
    self.assertEqual((0.,1.,0.), material.Ku1_atom_axis)
    self.assertEqual((1./sqrt(2.),1./sqrt(2.),0.), material.D_axis)
    self.assertEqual((1./sqrt(3.),1./sqrt(3.),1./sqrt(3.)), material.D_atom_axis)

  def test_param_default_initialization(self):
    material=magnumaf.Material()
    self.assertEqual(4e-7*pi, magnumaf.Constants.mu0           )
    self.assertEqual(221276.14886372554, magnumaf.Constants.gamma         )
    self.assertEqual(0., material.p             )
    self.assertEqual(0., material.D             )
    self.assertEqual(0., material.J_atom        )
    self.assertEqual(0., material.D_atom        )
    self.assertEqual(0., material.Ku1_atom        )

if __name__ == '__main__':
  unittest.main()
