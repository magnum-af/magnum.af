import unittest
import arrayfire as af
import magnumaf
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

  mesh=magnumaf.Mesh(nx, ny, nz, x/nx, y/ny, z/nz)
  material=magnumaf.Material()
  state=magnumaf.State(mesh, Ms = 0, m = m_af)
  
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

  def test_material_setter_and_getter_for_ms(self):
    Ms=1000
    state=magnumaf.State(magnumaf.Mesh(0,0,0,0,0,0), Ms = Ms, m = self.m_af)
    self.assertEqual(Ms, state.Ms)

  def test_mesh_setter_and_getter(self):
    state=magnumaf.State(magnumaf.Mesh(0,0,0,0,0,0), Ms = 0, m = self.m_af)
    mesh = magnumaf.Mesh(1,2,3,4,5,6)
    state.mesh=mesh
    self.assertEqual(1, state.mesh.n0)
    self.assertEqual(2, state.mesh.n1)
    self.assertEqual(3, state.mesh.n2)
    self.assertEqual(4, state.mesh.dx)
    self.assertEqual(5, state.mesh.dy)
    self.assertEqual(6, state.mesh.dz)

if __name__ == '__main__':
  unittest.main()
