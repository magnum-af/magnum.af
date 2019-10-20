import unittest
import arrayfire as af
import magnumaf
from numpy import zeros

class ZeemanTest(unittest.TestCase):
  dx = 1.e-10
  def test_zeeman_field(self):
    mesh=magnumaf.Mesh(3, 1, 1, self.dx, self.dx, self.dx)
    material=magnumaf.Material()
    m_np = zeros(3)
    m_np[0] = 1.0
    pystate=magnumaf.State(mesh, Ms = 0, m = af.reorder(af.from_ndarray(m_np), 1, 2, 3, 0))
    zeefield = af.constant(0.0, 3, 1, 1, 3, dtype=af.Dtype.f32)
    zeefield[0, 0, 0, 0]=1
    zeefield[1, 0, 0, 1]=1
    zeefield[2, 0, 0, 2]=1
    zee=magnumaf.ExternalField(zeefield)
    Llg=magnumaf.LLGIntegrator(0, [zee])

    af_heff = Llg.h(pystate)
    np_heff = af_heff.__array__()

    self.assertEqual(np_heff[0, 0, 0, 0], 1)
    self.assertEqual(np_heff[0, 0, 0, 1], 0)
    self.assertEqual(np_heff[0, 0, 0, 2], 0)
    self.assertEqual(np_heff[1, 0, 0, 0], 0)
    self.assertEqual(np_heff[1, 0, 0, 1], 1)
    self.assertEqual(np_heff[1, 0, 0, 2], 0)
    self.assertEqual(np_heff[2, 0, 0, 0], 0)
    self.assertEqual(np_heff[2, 0, 0, 1], 0)
    self.assertEqual(np_heff[2, 0, 0, 2], 1)

if __name__ == '__main__':
  unittest.main()
