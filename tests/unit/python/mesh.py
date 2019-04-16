import unittest
import arrayfire as af
from magnumaf import Mesh

class MeshTest(unittest.TestCase):
    n0 = 100
    n1 = 25
    n2 = 1
    dx = 5.e-7/100.
    dy = 1.25e-7/25.
    dz = 3.e-9

    def test_mesh_getter(self):
        mesh = Mesh(self.n0, self.n1, self.n2, self.dx, self.dy, self.dz)
        self.assertEqual(self.n0, mesh.n0)
        self.assertEqual(self.n1, mesh.n1)
        self.assertEqual(self.n2, mesh.n2)
        self.assertEqual(self.dx, mesh.dx)
        self.assertEqual(self.dy, mesh.dy)
        self.assertEqual(self.dz, mesh.dz)

    def test_mesh_setter_using_getter(self):
        n0 = 4*100
        n1 = 4*25
        n2 = 4*1
        dx = 4*5.e-7/100.
        dy = 4*1.25e-7/25.
        dz = 4*3.e-9
        mesh = Mesh(self.n0, self.n1, self.n2, self.dx, self.dy, self.dz)
        mesh.n0 = n0
        mesh.n1 = n1
        mesh.n2 = n2
        mesh.dx = dx
        mesh.dy = dy
        mesh.dz = dz
        self.assertEqual(n0, mesh.n0)
        self.assertEqual(n1, mesh.n1)
        self.assertEqual(n2, mesh.n2)
        self.assertEqual(dx, mesh.dx)
        self.assertEqual(dy, mesh.dy)
        self.assertEqual(dz, mesh.dz)

if __name__ == '__main__':
    unittest.main()
