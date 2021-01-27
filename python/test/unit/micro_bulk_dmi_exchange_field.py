import unittest
import arrayfire as af
from magnumaf import *

class BulkDMIExchangeFieldTest(unittest.TestCase):
    mesh = Mesh(4, 5, 6, 0.1, 0.2, 0.3)
    Ms = 8e5
    D = 1 # TODO better value
    A =1.3e-11
    def test_if_field_of_homog_m_is_zero(self):
        m = af.constant(0.0, self.mesh.nx, self.mesh.ny, self.mesh.nz, 3, dtype=af.Dtype.f64)
        m[:, :, :, 0] = 1
        state = State(self.mesh, self.Ms, m)
        dmi_exch = BulkDMIExchangeField(self.D, self.A)
        h = dmi_exch.H_in_Apm(state)
        h_wo_edges = h
        scalar_sum_over_h = af.sum(af.sum(af.sum(af.sum(h[1:-1, 1:-1, 1:-1], 3), 2), 1), 0).scalar()
        self.assertAlmostEqual(scalar_sum_over_h, 0.0, places = 12)

    # TODO check tilt at boundary vs analytical

if __name__ == '__main__':
    unittest.main()
