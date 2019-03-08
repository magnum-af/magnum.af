import arrayfire as af
from magnum_af import *
import sys
import time

af.info()

# Physical dimensions in [m]
x = 1e-9
y = 1e-9
z = 0.6e-9
# Discretization 
nx = 4
ny = 4
nz = 4

mesh = Mesh(nx, ny, nz, x/nx, y/ny, z/nz)
material = Material(ms = 1.4e6, A = 30e-12, alpha = 0.02, Ku1=1e4, Ku1_axis=[1,0,0])
m = Util.normed_homogeneous_field(nx, ny, nz, [1,0,0])
state = State(mesh, material, m)

polarization = Util.normed_homogeneous_field(nx, ny, nz, [0, 1, 0])

#demag = DemagField(mesh, material)
#exch = ExchangeField(mesh, material)
stt = SpinTransferTorqueField(polarization, nu_damp=.3, nu_field=.4, j_e=2e10)
aniso = UniaxialAnisotropyField(mesh, material)
Llg = LLGIntegrator(aniso)
#Llg = LLGIntegrator(stt, aniso)

print(stt.polarization_field)
print(state.m)

stream = open(sys.argv[1]+"m.dat", "w")
timer = time.time()
while state.t < 1e-9:
    Llg.llgstep(state)
    mean = af.mean(af.mean(af.mean(state.m, dim=0), dim=1), dim=2)
    print(state.t, mean[0,0,0,0].scalar(), mean[0,0,0,1].scalar(), mean[0,0,0,2].scalar())
    stream.write("%e, %e, %e, %e\n" %(state.t, mean[0,0,0,0].scalar(), mean[0,0,0,1].scalar(), mean[0,0,0,2].scalar()))
    stream.flush()
print("relaxed in", time.time() - timer, "[s]")
stream.close()


#heff = llg.get_fheff(state)
#print(heff[0,0,0,0].scalar())
#print(heff[0,0,0,1].scalar())
#print(heff[0,0,0,2].scalar())
