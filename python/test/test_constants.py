import unittest
import arrayfire as af
import magnumaf
from math import pi, sqrt


class ConstantsTest(unittest.TestCase):

    def test_constants_initialization(self):
        self.assertEqual(4e-7 * pi, magnumaf.Constants.mu0)
        self.assertEqual(221276.14886372554, magnumaf.Constants.gamma)


if __name__ == '__main__':
    unittest.main()
