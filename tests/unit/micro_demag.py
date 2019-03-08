import unittest
import arrayfire as af
import numpy as np
import magnum_af
import math

class MicroDemagTest(unittest.TestCase):
    nx=10
    dx=1.e-10
    def test_micro_demag_homogen_cube(self):
        mesh=magnum_af.Mesh(self.nx, self.nx, self.nx, self.dx, self.dx, self.dx)
        material=magnum_af.Material()
        material.ms = 1e5
        m=af.constant(0.0,self.nx, self.nx, self.nx, 3,dtype=af.Dtype.f64)
        m[:,:,:,0]=1.
        pystate=magnum_af.State(mesh,material,m)
        micro_demag=magnum_af.DemagField(mesh,material)
        Llg=magnum_af.LLGIntegrator([micro_demag])
        self.assertAlmostEqual(Llg.get_E(pystate), 1./6. * (self.nx*self.dx)**3 * material.ms**2 * magnum_af.Constants.mu0)

if __name__ == '__main__':
  unittest.main()
