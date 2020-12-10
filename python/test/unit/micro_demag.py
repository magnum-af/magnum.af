import unittest
import arrayfire as af
import numpy as np
import magnumaf
import math

class MicroDemagTest(unittest.TestCase):
    nx=10
    dx=1.e-10
    def test_micro_demag_homogen_cube(self):
        mesh=magnumaf.Mesh(self.nx, self.nx, self.nx, self.dx, self.dx, self.dx)
        m=af.constant(0.0, self.nx, self.nx, self.nx, 3, dtype=af.Dtype.f64)
        m[:, :, :, 0]=1.
        state=magnumaf.State(mesh, Ms = 1e5, m = m)
        micro_demag=magnumaf.DemagField(mesh)
        llg=magnumaf.LLGIntegrator(alpha = 0, terms = [micro_demag])
        self.assertAlmostEqual(llg.E(state), 1./6. * (self.nx*self.dx)**3 * state.Ms**2 * magnumaf.Constants.mu0)

if __name__ == '__main__':
  unittest.main()
