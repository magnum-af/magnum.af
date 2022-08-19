import unittest
import arrayfire as af
from magnumaf import *

# Note: this is just a wrapping test/demo, numeric tests are in c++


class MicroExchangePBCTest(unittest.TestCase):
    def test_homogen_m(self):
        mesh = Mesh(8, 10, 12, 1e-9, 2e-9, 3e-9)
        m = af.constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, dtype=af.Dtype.f64)
        m[:, :, :, 2] = 1.
        state = State(mesh, Ms=1e5, m=m)
        A = 1e-12
        exch_pbc = ExchangeFieldPBC(A, mesh)
        H = exch_pbc.H_in_Apm(state)
        H_sum = af.sum(af.sum(af.sum(af.sum(H, 0), 1), 2), 3).scalar()
        self.assertAlmostEqual(H_sum, 0, 5)


if __name__ == '__main__':
    unittest.main()
