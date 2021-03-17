import arrayfire as af
from magnumaf import *

# Physical dimensions in [m]
x = 1e-9
y = 1e-9
z = 1e-9
# Discretization
nx = 1
ny = 1
nz = 1

polarization = Magnetization.homogeneous(nx, ny, nz, [0, 1, 0])
Util.test_sum_of_difference_of_abs(polarization, polarization)

fieldlike = SpinTransferTorqueField(polarization, 1., 1., 1.)
Util.test_sum_of_difference_of_abs(polarization, fieldlike.polarization_field)
fieldlike.polarization_field = af.constant(2.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
Util.test_sum_of_difference_of_abs(2*polarization, fieldlike.polarization_field)
fieldlike.polarization_field = polarization

m = Magnetization.homogeneous(nx, ny, nz, [1, 0, 0])
print(fieldlike.polarization_field)
print(m)

state = State(Mesh(nx, ny, nz, x/nx, y/ny, z/nz), Material(ms = 8e5, A = 1.3e-11, alpha = 1.), m)
llg = LLGIntegrator([fieldlike])
heff = llg.h(state)
print(heff[0, 0, 0, 0].scalar())
print(heff[0, 0, 0, 1].scalar())
print(heff[0, 0, 0, 2].scalar())
