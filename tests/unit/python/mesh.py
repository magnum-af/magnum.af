import unittest
import arrayfire as af
from magnumaf import Mesh

class MeshTest(unittest.TestCase):
    nx = 100
    ny = 25
    nz = 1
    dx = 5.e-7/100.
    dy = 1.25e-7/25.
    dz = 3.e-9

    def test_mesh_getter(self):
        mesh = Mesh(self.nx, self.ny, self.nz, self.dx, self.dy, self.dz)
        self.assertEqual(self.nx, mesh.nx)
        self.assertEqual(self.ny, mesh.ny)
        self.assertEqual(self.nz, mesh.nz)
        self.assertEqual(self.dx, mesh.dx)
        self.assertEqual(self.dy, mesh.dy)
        self.assertEqual(self.dz, mesh.dz)

    def test_mesh_setter_using_getter(self):
        nx = 4*100
        ny = 4*25
        nz = 4*1
        dx = 4*5.e-7/100.
        dy = 4*1.25e-7/25.
        dz = 4*3.e-9
        mesh = Mesh(self.nx, self.ny, self.nz, self.dx, self.dy, self.dz)
        mesh.nx = nx
        mesh.ny = ny
        mesh.nz = nz
        mesh.dx = dx
        mesh.dy = dy
        mesh.dz = dz
        self.assertEqual(nx, mesh.nx)
        self.assertEqual(ny, mesh.ny)
        self.assertEqual(nz, mesh.nz)
        self.assertEqual(dx, mesh.dx)
        self.assertEqual(dy, mesh.dy)
        self.assertEqual(dz, mesh.dz)

if __name__ == '__main__':
    unittest.main()
