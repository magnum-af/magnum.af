# run with:
#PYTHONPATH=../../build/python/ python sp4_output.py path/to/write/in
import arrayfire as af
import magnumaf
import os
import sys
import matplotlib.pyplot as plt

path = sys.argv[1]
if os.path.exists(path) is False:
    print "path does not exist!"
    exit()
af.info()

meshvar=magnumaf.Mesh(  100, 25, 1, 5.e-7/100, 1.25e-7/25, 3.e-9)
m=af.constant(0.0, 100, 25, 1, 3, dtype=af.Dtype.f64)

material=magnumaf.Material()
state.Ms    (8e5)
material.A     (1.3e-11)
material.alpha (1)

m[1:-1, :, :, 0] = af.constant(1.0, 100-2, 25, 1, 1, dtype=af.Dtype.f64);
m[0, :, :, 1]    = af.constant(1.0, 1    , 25, 1, 1, dtype=af.Dtype.f64);
m[-1, :, :, 1]   = af.constant(1.0, 1    , 25, 1, 1, dtype=af.Dtype.f64);
pystate=magnumaf.State(meshvar, material, m)
pystate.write_vti(path+"minit")

demag=magnumaf.DemagField(meshvar, material)
exch=magnumaf.ExchangeField(meshvar, material)
llg=magnumaf.LLGIntegrator([pystate, demag, exch])

print "relax --------------------"
while pystate.t() < 1e-9:
  llg.step(pystate)
pystate.write_vti(path+"mrelax")

print "switch --------------------"
llg.set_state0_alpha(0.02)# this should be changed in cpp version

zeeswitch = af.constant(0.0, 1, 1, 1, 3, dtype=af.Dtype.f64)
zeeswitch[0, 0, 0, 0]=-24.6e-3/material.print_mu0()
zeeswitch[0, 0, 0, 1]=+4.3e-3/material.print_mu0()
zeeswitch[0, 0, 0, 2]=0.0
zeeswitch = af.tile(zeeswitch, 100, 25, 1)
zee=magnumaf.ExternalField(zeeswitch)
llg.add_terms(zee)

with open(path + 'm.dat', 'w') as f:
  while pystate.t() < 2e-9:
    llg.step(pystate)
    f.write("%10.12f %10.12f %10.12f %10.12f\n" % (pystate.t(), pystate.mean_m(0), pystate.mean_m(1), pystate.mean_m(2)))

# pyplot fails:
#with open(path+'m.dat') as f:
#  lines = f.readlines()
#  t = [line.split()[0] for line in lines]
#  x = [line.split()[1] for line in lines]
#  y = [line.split()[2] for line in lines]
#  z = [line.split()[3] for line in lines]
#plt.plot(t, y)
#plt.show()
