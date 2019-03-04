import arrayfire as af
from magnum_af import *

# Discretization 
nx = 1
ny = 1
nz = 1

def test(a, b, verbose = True):
    c = af.sum(af.sum(af.sum(af.sum(af.abs(a)-af.abs(b),0),1),2),3).scalar()
    if (c != 0.):
        if (verbose == True):
            print ("Error")
        return False
    else:
        if (verbose == True):
            print ("Success")
        return True

polarization = af.constant(1.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
test(polarization, polarization)

fieldlike = SpinTransferTorqueField(polarization, 1., 1., 1.)
test(polarization, fieldlike.polarization_field)
fieldlike.polarization_field = af.constant(2.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
test(2*polarization, fieldlike.polarization_field)
