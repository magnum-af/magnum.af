# run with: 
#PYTHONPATH=../../build/src/ python sp4_output.py path/to/write/in
import arrayfire as af
import magnum_af
import os
import sys
import matplotlib.pyplot as plt

path = sys.argv[1]
if os.path.exists(path) is False:
    print "path does not exist!"
    exit()
af.info()

meshvar=magnum_af.Mesh(  100,25,1,5.e-7/100,1.25e-7/25,3.e-9)
m=af.constant(0.0,100,25,1,3,dtype=af.Dtype.f64)

material=magnum_af.Material()
material.ms    (8e5)
material.A     (1.3e-11)
material.alpha (1)

m[1:-1,:,:,0] = af.constant(1.0,100-2,25,1,1,dtype=af.Dtype.f64);
m[0,:,:,1]    = af.constant(1.0,1    ,25,1,1,dtype=af.Dtype.f64);
m[-1,:,:,1]   = af.constant(1.0,1    ,25,1,1,dtype=af.Dtype.f64);
pystate=magnum_af.State(meshvar,material,m)
pystate.py_vti_writer_micro(path+"minit")

demag=magnum_af.DemagField(meshvar,material)
exch=magnum_af.ExchangeField(meshvar,material)
Llg=magnum_af.LLGIntegrator([pystate,demag,exch])

print "relax --------------------"
while pystate.t() < 1e-9:
  Llg.step(pystate)
pystate.py_vti_writer_micro(path+"mrelax")

print "switch --------------------"
Llg.set_state0_alpha(0.02)# this should be changed in cpp version

zeeswitch = af.constant(0.0,1,1,1,3,dtype=af.Dtype.f64)
zeeswitch[0,0,0,0]=-24.6e-3/material.print_mu0()
zeeswitch[0,0,0,1]=+4.3e-3/material.print_mu0()
zeeswitch[0,0,0,2]=0.0
zeeswitch = af.tile(zeeswitch,100,25,1)
zee=magnum_af.ExternalField(zeeswitch)
Llg.add_terms(zee)

with open(path + 'm.dat', 'w') as f:
  while pystate.t() < 2e-9:
    Llg.step(pystate)
    f.write("%10.12f %10.12f %10.12f %10.12f\n" % (pystate.t(), pystate.meanxyz(0), pystate.meanxyz(1), pystate.meanxyz(2)))

# pyplot fails:
#with open(path+'m.dat') as f:
#  lines = f.readlines()
#  t = [line.split()[0] for line in lines]
#  x = [line.split()[1] for line in lines]
#  y = [line.split()[2] for line in lines]
#  z = [line.split()[3] for line in lines]
#plt.plot(t,y)
#plt.show()
