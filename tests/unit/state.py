import unittest
import arrayfire as af
import magnum_af
from numpy import zeros, ones
from math import pi, sqrt

class StateTest(unittest.TestCase):

  x = 1.e-9
  y = 1.e-9
  z = 1.e-9
  nx = 1
  ny = 1
  nz = 1

  m_np = zeros(3)
  m_np[0] = 1.0
  m_af = af.reorder(af.from_ndarray(m_np),1,2,3,0)

  mesh=magnum_af.Mesh(nx, ny, nz, x/nx, y/ny, z/nz)
  param=magnum_af.Param()
  state=magnum_af.State(mesh, param, m_af)
  
  def test_state_initialization(self):
    self.assertEqual(1., self.state.m[0,0,0,0].scalar())
    self.assertEqual(0., self.state.m[0,0,0,1].scalar())
    self.assertEqual(0., self.state.m[0,0,0,2].scalar())

  def test_state_set_m(self):
    m_set_np = zeros(3)
    m_set_np[1] = 1.0
    m_set_af = af.from_ndarray(m_set_np)
    m_set_af = af.reorder(m_set_af,1,2,3,0)

    self.state.m = m_set_af
    self.assertEqual(0., self.state.m[0,0,0,0].scalar())
    self.assertEqual(1., self.state.m[0,0,0,1].scalar())
    self.assertEqual(0., self.state.m[0,0,0,2].scalar())

  def test_state_set_wrong_m(self): 
    m_np = ones(3)
    m_af = af.from_ndarray(m_np)
    m_af = af.reorder(m_af,1,2,3,0)

    print("The following warning is provoked:")
    self.state.m = m_af
    self.assertEqual(1., self.state.m[0,0,0,0].scalar())
    self.assertEqual(1., self.state.m[0,0,0,1].scalar())
    self.assertEqual(1., self.state.m[0,0,0,2].scalar())

if __name__ == '__main__':
  unittest.main()
