import unittest
import arrayfire as af
from magnumaf import CubicAnisotropyField, State, Mesh, Constants

class StateTest(unittest.TestCase):
    def test_Kc1_energy(self):
        Kc1 = 1e6
        caniso = CubicAnisotropyField(Kc1, 0, 0, [1, 0, 0], [0, 1, 0])
        mesh = Mesh(1, 1, 1, 1e-9, 1e-9, 1e-9)
        m = af.constant(0, 2, 1, 1, 3, dtype = af.Dtype.f64)
        m[0, :, :, 0] = 1.
        m[0, :, :, 1] = 1.
        state = State(mesh, Ms = 1., m = m)
        E = caniso.E(state)
        # 250000 from: Kc1_E_density_analytic = Kc1 * (c1m2 * c2m2 + c1m2 * c3m2 + c2m2 * c3m2);
        self.assertAlmostEqual(E, 250000 * mesh.dx * mesh.dy * mesh.dz, 36)

if __name__ == '__main__':
  unittest.main()
