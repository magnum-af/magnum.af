import unittest
import arrayfire as af
import numpy as np
import magnumaf
import math


class AtomisticDmiFieldTest(unittest.TestCase):
    # arbitrary parameters:
    p = 1e-20
    D = 1e20
    dx = 1.1
    # TODO init test
    # TODO#D_atom = 8.5
    # TODO#D_atom_axis = [1., 1., 1.]
    # TODO#self.assertEqual(self.D_atom, D_atom)
    # TODO#self.assertEqual((1./sqrt(3.), 1./sqrt(3.), 1./sqrt(3.)), D_atom_axis)

    def test_atomistic_dmi_2_1_1_z_z(self):
        mesh = magnumaf.Mesh(2, 1, 1, self.dx, self.dx, self.dx)
        p = self.p
        D_atom = self.D
        m = af.constant(0.0, 2, 1, 1, 3, dtype=af.Dtype.f64)
        m[0, 0, 0, 0] = 0
        m[0, 0, 0, 1] = 0
        m[0, 0, 0, 2] = 1

        m[1, 0, 0, 0] = 0
        m[1, 0, 0, 1] = 0
        m[1, 0, 0, 2] = 1

        state = magnumaf.State(mesh, Ms=p, m=m)
        atom_dmi = magnumaf.AtomisticDmiField(D_atom, [0, 0, 1])

        self.assertAlmostEqual(atom_dmi.Energy_in_J(state), 0)

        af_heff = atom_dmi.H_in_Apm(state)
        np_heff = af_heff.__array__()

        self.assertAlmostEqual(
            np_heff[0, 0, 0, 0], -D_atom/magnumaf.Constants.mu0/p)
        self.assertLess(math.fabs(
            (np_heff[0, 0, 0, 0] - (-D_atom/magnumaf.Constants.mu0/p))/np_heff[0, 0, 0, 0]), 1e-15)
        self.assertAlmostEqual(np_heff[0, 0, 0, 1], 0)
        self.assertAlmostEqual(np_heff[0, 0, 0, 2], 0)
        self.assertAlmostEqual(
            np_heff[1, 0, 0, 0], D_atom/magnumaf.Constants.mu0/p)
        self.assertLess(math.fabs(
            (np_heff[1, 0, 0, 0] - D_atom/magnumaf.Constants.mu0/p)/np_heff[1, 0, 0, 0]), 1e-15)
        self.assertAlmostEqual(np_heff[1, 0, 0, 1], 0)
        self.assertAlmostEqual(np_heff[1, 0, 0, 2], 0)

    def test_atomistic_dmi_2_1_1_z_x(self):
        mesh = magnumaf.Mesh(2, 1, 1, self.dx, self.dx, self.dx)
        p = self.p
        D_atom = self.D
        m = af.constant(0.0, 2, 1, 1, 3, dtype=af.Dtype.f64)
        m[0, 0, 0, 0] = 0
        m[0, 0, 0, 1] = 0
        m[0, 0, 0, 2] = 1

        m[1, 0, 0, 0] = 1
        m[1, 0, 0, 1] = 0
        m[1, 0, 0, 2] = 0

        state = magnumaf.State(mesh, Ms=p, m=m)
        atom_dmi = magnumaf.AtomisticDmiField(D_atom, [0, 0, 1])

        self.assertAlmostEqual(atom_dmi.Energy_in_J(state), - D_atom)

        af_heff = atom_dmi.H_in_Apm(state)
        np_heff = af_heff.__array__()

        self.assertAlmostEqual(np_heff[0, 0, 0, 0], 0)
        self.assertAlmostEqual(np_heff[0, 0, 0, 1], 0)
        self.assertAlmostEqual(
            np_heff[0, 0, 0, 2], D_atom/magnumaf.Constants.mu0/p)
        self.assertLess(math.fabs(
            (np_heff[0, 0, 0, 2] - D_atom/magnumaf.Constants.mu0/p)/np_heff[0, 0, 0, 2]), 1e-15)

        self.assertAlmostEqual(
            np_heff[1, 0, 0, 0], D_atom/magnumaf.Constants.mu0/p)
        self.assertLess(math.fabs(
            (np_heff[1, 0, 0, 0] - D_atom/magnumaf.Constants.mu0/p)/np_heff[1, 0, 0, 0]), 1e-15)
        self.assertAlmostEqual(np_heff[1, 0, 0, 1], 0)
        self.assertAlmostEqual(np_heff[1, 0, 0, 2], 0)

    def test_atomistic_dmi_1_1_2_z_z(self):
        mesh = magnumaf.Mesh(1, 1, 2, self.dx, self.dx, self.dx)
        p = self.p
        D_atom = self.D
        m = af.constant(0.0, 1, 1, 2, 3, dtype=af.Dtype.f64)
        m[0, 0, 0, 0] = 0
        m[0, 0, 0, 1] = 0
        m[0, 0, 0, 2] = 1

        m[0, 0, 1, 0] = 1
        m[0, 0, 1, 1] = 0
        m[0, 0, 1, 2] = 0

        state = magnumaf.State(mesh, Ms=p, m=m)
        atom_dmi = magnumaf.AtomisticDmiField(D_atom, [0, 0, 1])

        self.assertAlmostEqual(atom_dmi.Energy_in_J(state), 0)

        af_heff = atom_dmi.H_in_Apm(state)
        np_heff = af_heff.__array__()

        self.assertAlmostEqual(np_heff[0, 0, 0, 0], 0)
        self.assertAlmostEqual(np_heff[0, 0, 0, 1], 0)
        self.assertAlmostEqual(np_heff[0, 0, 0, 2], 0)
        self.assertAlmostEqual(np_heff[0, 0, 1, 0], 0)
        self.assertAlmostEqual(np_heff[0, 0, 1, 1], 0)
        self.assertAlmostEqual(np_heff[0, 0, 1, 2], 0)


if __name__ == '__main__':
    unittest.main()
