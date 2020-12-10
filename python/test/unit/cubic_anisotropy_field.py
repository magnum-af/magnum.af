import unittest
import arrayfire as af
from magnumaf import CubicAnisotropyField, State, Mesh, Constants

class StateTest(unittest.TestCase):
    def test_Kc1_energy(self):
        Kc1 = 1e6
        caniso = CubicAnisotropyField(Kc1, 0, 0, [1, 0, 0], [0, 1, 0])
        mesh = Mesh(1, 1, 1, 1e-9, 1e-9, 1e-9)
        m = af.constant(0, mesh.nx, mesh.ny, mesh.nz, 3, dtype = af.Dtype.f64)
        m[:, :, :, 0] = 1.
        m[:, :, :, 1] = 1.
        state = State(mesh, Ms = 1., m = m)
        E = caniso.Energy_in_J(state)
        # 250000 from: Kc1_E_density_analytic = Kc1 * (c1m2 * c2m2 + c1m2 * c3m2 + c2m2 * c3m2);
        self.assertAlmostEqual(E, 250000 * mesh.dx * mesh.dy * mesh.dz, 36)

    def test_Kc1_energy_array_input(self):
        mesh = Mesh(2, 2, 2, 1e-9, 1e-9, 1e-9)
        Kc1_val = 1e6
        Kc1 = af.constant(Kc1_val, mesh.nx, mesh.ny, mesh.nz, 1, dtype = af.Dtype.f64)
        Kc2 = af.constant(0, mesh.nx, mesh.ny, mesh.nz, 1, dtype = af.Dtype.f64)
        Kc3 = af.constant(0, mesh.nx, mesh.ny, mesh.nz, 1, dtype = af.Dtype.f64)
        c1 = af.constant(0, mesh.nx, mesh.ny, mesh.nz, 3, dtype = af.Dtype.f64)
        c1[:, :, :, 0] = 1
        c2 = af.constant(0, mesh.nx, mesh.ny, mesh.nz, 3, dtype = af.Dtype.f64)
        c2[:, :, :, 1] = 1
        caniso = CubicAnisotropyField(Kc1, Kc2, Kc3, c1, c2)
        m = af.constant(0, mesh.nx, mesh.ny, mesh.nz, 3, dtype = af.Dtype.f64)
        m[:, :, :, 0] = 1.
        m[:, :, :, 1] = 1.
        state = State(mesh, Ms = 1., m = m)
        E = caniso.Energy_in_J(state)
        # 250000 from: Kc1_E_density_analytic = Kc1 * (c1m2 * c2m2 + c1m2 * c3m2 + c2m2 * c3m2);
        self.assertAlmostEqual(E, mesh.nx * mesh.ny * mesh.nz * 250000 * mesh.dx * mesh.dy * mesh.dz, 20) # precision 36 for nx=ny=nz=1

if __name__ == '__main__':
  unittest.main()
