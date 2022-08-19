import unittest
import arrayfire as af
import magnumaf


class FunctionsTest(unittest.TestCase):

    def test_spacial_mean_in_region(self):
        vectorfield = af.constant(0, 6, 1, 1, 3, dtype=af.Dtype.f64)
        vectorfield[0:3, 0, 0, 0] = 2.5
        vectorfield[0:3, 0, 0, 1] = 3.5
        vectorfield[0:3, 0, 0, 2] = 4.5
        region = af.constant(0, 6, 1, 1, 1, dtype=af.Dtype.f64)
        region[0:3] = 1
        res = magnumaf.Util.spacial_mean_in_region(vectorfield, region)
        self.assertEqual(res[0], 2.5)
        self.assertEqual(res[1], 3.5)
        self.assertEqual(res[2], 4.5)


if __name__ == '__main__':
    unittest.main()
