import unittest
import arrayfire as af
import magnumaf
from numpy import zeros, ones
from math import pi, sqrt

verbose = False
# verbose = True  # enable for dev

class StateTest(unittest.TestCase):

  x = 1.e-9
  y = 1.e-9
  z = 1.e-9
  nx = 1
  ny = 1
  nz = 1

  m_np = zeros(3)
  m_np[0] = 1.0
  m_af = af.reorder(af.from_ndarray(m_np), 1, 2, 3, 0)

  mesh=magnumaf.Mesh(nx, ny, nz, x/nx, y/ny, z/nz)
  state=magnumaf.State(mesh, Ms = 0, m = m_af)

  def test_state_initialization(self):
    self.assertEqual(1., self.state.m[0, 0, 0, 0].scalar())
    self.assertEqual(0., self.state.m[0, 0, 0, 1].scalar())
    self.assertEqual(0., self.state.m[0, 0, 0, 2].scalar())

  def test_state_set_m(self):
    m_set_np = zeros(3)
    m_set_np[1] = 1.0
    m_set_af = af.from_ndarray(m_set_np)
    m_set_af = af.reorder(m_set_af, 1, 2, 3, 0)

    self.state.m = m_set_af
    self.assertEqual(0., self.state.m[0, 0, 0, 0].scalar())
    self.assertEqual(1., self.state.m[0, 0, 0, 1].scalar())
    self.assertEqual(0., self.state.m[0, 0, 0, 2].scalar())

  def test_state_normalize_m(self):
    m_np = ones(3)
    m_af = af.from_ndarray(m_np)
    m_af = af.reorder(m_af, 1, 2, 3, 0)

    self.state.m = m_af
    self.assertEqual(1./sqrt(3), self.state.m[0, 0, 0, 0].scalar())
    self.assertEqual(1./sqrt(3), self.state.m[0, 0, 0, 1].scalar())
    self.assertEqual(1./sqrt(3), self.state.m[0, 0, 0, 2].scalar())

  def test_material_setter_and_getter_for_ms(self):
    Ms=1000
    state=magnumaf.State(magnumaf.Mesh(0, 0, 0, 0, 0, 0), Ms = Ms, m = self.m_af)
    self.assertEqual(Ms, state.Ms)

  def test_material_setter_and_getter_for_Ms_field(self):
    Msval = 1001.0
    Ms_field = af.constant(Msval, 1, 1, 1, 1, af.Dtype.f64)
    state=magnumaf.State(self.mesh, Ms = Ms_field, m = self.m_af)
    if verbose: print(state.Ms)
    self.assertEqual(Msval, state.Ms.scalar()) # test getter
    state.Ms = af.constant(Msval + 1.5, 1, 1, 1, 1, af.Dtype.f64) # setter
    if verbose: print(state.Ms)
    self.assertEqual(Msval + 1.5, state.Ms.scalar()) # test setter-getter

  # # Note: Ms switching between float and af.Array is not supported
  # def test_material_setter_and_getter_for_Ms_field_scalar_init(self):
  #   Msval = 1001
  #   Ms_field = af.constant(Msval, 1, 1, 1, 1, af.Dtype.f64)
  #   state=magnumaf.State(self.mesh, Ms = 0.0, m = self.m_af)
  #   if verbose: print(state.Ms)
  #   state.Ms = Ms_field # setter
  #   if verbose: print(state.Ms)
  #   self.assertEqual(Msval, state.Ms.scalar())

  def test_mesh_setter_and_getter(self):
    state=magnumaf.State(magnumaf.Mesh(0, 0, 0, 0, 0, 0), Ms = 0, m = self.m_af)
    mesh = magnumaf.Mesh(1, 2, 3, 4, 5, 6)
    state.mesh=mesh
    self.assertEqual(1, state.mesh.nx)
    self.assertEqual(2, state.mesh.ny)
    self.assertEqual(3, state.mesh.nz)
    self.assertEqual(4, state.mesh.dx)
    self.assertEqual(5, state.mesh.dy)
    self.assertEqual(6, state.mesh.dz)
  def test_mean_M(self):
    Ms = 1e8
    Ms_field = af.constant(Ms, 1, 1, 1, 1)
    m0 = af.constant(0, 1, 1, 1, 3)
    m0[0, 0, 0, 0] = 1 # x dir
    state=magnumaf.State(self.mesh, Ms = Ms_field, m = m0)
    Mx, My, Mz = state.mean_M()
    self.assertEqual(Ms, Mx)
    self.assertEqual(0, My)
    self.assertEqual(0, Mz)
  def test_mean_M_no_Ms_field(self):
    Ms = 1e8
    m0 = af.constant(0, 1, 1, 1, 3)
    m0[0, 0, 0, 0] = 1 # x dir
    state=magnumaf.State(self.mesh, Ms = Ms, m = m0)
    Mx, My, Mz = state.mean_M()
    self.assertEqual(Ms, Mx)
    self.assertEqual(0, My)
    self.assertEqual(0, Mz)


class StateMsDimTest(unittest.TestCase):
  m = af.constant(0, 1, 1, 1, 3, af.Dtype.f64)
  m[:, :, :, 0] = 1.

  Ms = af.constant(1, 1, 1, 1, 1, af.Dtype.f64)
  state=magnumaf.State(magnumaf.Mesh(1, 1, 1, 1, 1, 1), Ms, m)

  def test_Ms_dim(self):
      Ms_field = self.state.Ms_field
      self.assertEqual(Ms_field.dims(), self.Ms.dims())


class StateMsLegacyDimTest(unittest.TestCase):
  m = af.constant(0, 1, 1, 1, 3, af.Dtype.f64)
  m[:, :, :, 0] = 1.

  Ms_legacy = af.constant(1, 1, 1, 1, 3, af.Dtype.f64)
  state_legacy=magnumaf.State(magnumaf.Mesh(1, 1, 1, 1, 1, 1), Ms_legacy, m)

  def test_legacy_Ms_dim(self):
      Ms_field = self.state_legacy.Ms_field
      self.assertNotEqual(Ms_field.dims(), self.Ms_legacy.dims())
      self.assertEqual(Ms_field.dims()[0], 1)
      self.assertEqual(len(Ms_field.dims()), 1)


class State_mean_m_Test(unittest.TestCase):
  """2x2x2 Cube where m0 points in x-dir and is half empty. Zero cells should be neglected in mean_m calculation"""
  expected_mean_m = [1, 0, 0]
  Ms = 8e5
  mesh=magnumaf.Mesh(2, 2, 2, 1e-9, 1e-9, 1e-9)
  m0 = af.constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, af.Dtype.f64)
  m0[0, :, :, 0] = 1
  # print(m0)

  def test_scalar_Ms(self):
    state=magnumaf.State(self.mesh, Ms = self.Ms, m = self.m0)
    mx, my, mz = state.mean_m()
    self.assertEqual(mx, self.expected_mean_m[0])
    self.assertEqual(my, self.expected_mean_m[1])
    self.assertEqual(mz, self.expected_mean_m[2])

  def test_cell_wise_Ms(self):
    Ms_field = self.Ms * self.m0[:, :, :, 0]
    state=magnumaf.State(self.mesh, Ms = Ms_field, m = self.m0)
    mx, my, mz = state.mean_m()
    self.assertEqual(mx, self.expected_mean_m[0])
    self.assertEqual(my, self.expected_mean_m[1])
    self.assertEqual(mz, self.expected_mean_m[2])


if __name__ == '__main__':
  unittest.main()
