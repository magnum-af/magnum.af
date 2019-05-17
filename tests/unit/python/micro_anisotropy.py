import unittest
import arrayfire as af
from math import sqrt
from magnumaf import UniaxialAnisotropyField

class MicroAnisotropyTest(unittest.TestCase):
    def test_Ku1_axis_getter_2_0_0(self):
        Ku1 = 2
        Ku1_axis = [2, 0, 0]
        aniso = UniaxialAnisotropyField(Ku1, Ku1_axis)
        Ku1_axis_renormed = [x/(sqrt(Ku1_axis[0]**2 + Ku1_axis[1]**2 + Ku1_axis[2]**2)) for x in Ku1_axis]
        self.assertEqual(aniso.Ku1, Ku1)
        self.assertEqual(aniso.Ku1_axis[0], Ku1_axis_renormed[0])
        self.assertEqual(aniso.Ku1_axis[1], Ku1_axis_renormed[1])
        self.assertEqual(aniso.Ku1_axis[2], Ku1_axis_renormed[2])

    def test_Ku1_axis_getter_1_1_1(self):
        Ku1 = 1
        Ku1_axis = [1, 1, 1]
        aniso = UniaxialAnisotropyField(Ku1, Ku1_axis)
        Ku1_axis_renormed = [x/(sqrt(Ku1_axis[0]**2 + Ku1_axis[1]**2 + Ku1_axis[2]**2)) for x in Ku1_axis]
        self.assertEqual(aniso.Ku1, Ku1)
        self.assertEqual(aniso.Ku1_axis[0], Ku1_axis_renormed[0])
        self.assertEqual(aniso.Ku1_axis[1], Ku1_axis_renormed[1])
        self.assertEqual(aniso.Ku1_axis[2], Ku1_axis_renormed[2])

    def test_Ku1_axis_getter_default(self):
        Ku1 = 1
        default = [0, 0, 1]
        aniso = UniaxialAnisotropyField(Ku1)
        self.assertEqual(aniso.Ku1, Ku1)
        self.assertEqual(aniso.Ku1_axis[0], default[0])
        self.assertEqual(aniso.Ku1_axis[1], default[1])
        self.assertEqual(aniso.Ku1_axis[2], default[2])

if __name__ == '__main__':
    unittest.main()
