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
  mode = 6

  Ku1_axis = [2.0, 0., 0.]
  Ku1_atom_axis = [0., 3., 0.]
  D_axis = [1., 1., 0.]
  D_atom_axis = [1., 1., 1.]

  def test_param_initialization(self):
    param=magnum_af.pyParam(alpha=self.alpha, T=self.T, ms=self.ms, A=self.A, D=self.D, Ku1=self.Ku1, J_atom=self.J_atom, D_atom=self.D_atom, mode=self.mode, Ku1_axis=self.Ku1_axis, Ku1_atom_axis=self.Ku1_atom_axis, D_axis=self.D_axis, D_atom_axis=self.D_atom_axis)

    self.assertEqual(self.alpha, param.get_alpha())
    self.assertEqual(self.T, param.get_T())
    self.assertEqual(self.ms, param.get_ms())
    self.assertEqual(self.A, param.get_A())
    self.assertEqual(self.D, param.get_D())
    self.assertEqual(self.Ku1, param.get_Ku1())
    self.assertEqual(self.J_atom, param.get_J_atom())
    self.assertEqual(self.D_atom, param.get_D_atom())
    self.assertEqual(self.mode, param.get_mode())
    self.assertEqual((1.,0.,0.), param.get_Ku1_axis())
    self.assertEqual((0.,1.,0.), param.get_Ku1_atom_axis())
    self.assertEqual((1./sqrt(2.),1./sqrt(2.),0.), param.get_D_axis())
    self.assertEqual((1./sqrt(3.),1./sqrt(3.),1./sqrt(3.)), param.get_D_atom_axis())

  def test_param_default_initialization(self):
    param=magnum_af.pyParam()
    self.assertEqual(4e-7*pi, param.get_mu0()           )
    self.assertEqual(221276.1488637255, param.get_gamma()         )
    self.assertEqual((0.,0.,1.), param.get_Ku1_axis())
    self.assertEqual(0., param.get_Ku1_axis_x()    )
    self.assertEqual(0., param.get_Ku1_axis_y()    )
    self.assertEqual(1., param.get_Ku1_axis_z()    )
    self.assertEqual(0., param.get_Ku1_atom_axis_x() )
    self.assertEqual(0., param.get_Ku1_atom_axis_y() )
    self.assertEqual(1., param.get_Ku1_atom_axis_z() )
    self.assertEqual(0., param.get_ms()            )
    self.assertEqual(0., param.get_alpha()         )
    self.assertEqual(0., param.get_A()             )
    self.assertEqual(0., param.get_p()             )
    self.assertEqual(6 , param.get_mode()          )
    self.assertEqual(0., param.get_D()             )
    self.assertEqual(0., param.get_Ku1()           )
    self.assertEqual(0., param.get_J_atom()        )
    self.assertEqual(0., param.get_D_atom()        )
    self.assertEqual(0., param.get_Ku1_atom()        )

if __name__ == '__main__':
  unittest.main()
