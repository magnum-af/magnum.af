#!/usr/bin/python3
# Credit: Florian Bruckner, Amil Ducevic
# [1] F. Bruckner et al., “Strayfield calculation for micromagnetic simulations using true periodic boundary conditions,” Scientific Reports, vol. 11, p. 9202, Dec. 2021.
# expects N as first positional argument
# PBC demag is replaced by regular demag if secod pos argumet is "NOPBC"

import arrayfire as af
import numpy as np
from magnumaf import *
import time

args = parse()
N = int(args.posargs[0])
use_true_pbc = False if args.posargs[1] == "NOPBC" else True
print("use_true_pbc =", use_true_pbc, flush = True)

f = 1e8

k_dir = np.array(((1.688281290426342229e-01, -7.18775080428410873e-01, 6.744326850020654351e-01),
                  (-1.24488794472133734e-02, 8.827421556265807601e-01, 4.696927847862328309e-01),
                  (4.360199027887372986e-02, 9.960452028850761419e-01, -7.74133079860866157e-02),
                  (3.749512630383480816e-01, 6.881354883832631053e-01, -6.21193287128471327e-01),
                  (8.860065554198717219e-01, 2.300892602270045995e-01, -4.02555978816864223e-01),
                  (-7.89116308031293234e-02, -7.44290764462179010e-01, -6.63177361239280949e-01),
                  (1.534426618908853179e-01, 6.661409155280296757e-01, 7.298709681658244186e-01),
                  (3.629777351625293469e-02, 7.384491837339980380e-01, 6.733314745950575997e-01),
                  (-6.22338299750334056e-01, -4.54495719851937984e-01, -6.37282261874698607e-01),
                  (5.091810893444697061e-01, -7.86185846048253700e-01, 3.502091285608797122e-01),
                  (-7.42309465286212555e-01, 5.245654409560085440e-01, 4.169025736321682052e-01),
                  (2.976579786344664136e-01, 5.698206235530806074e-01, -7.65966177274703174e-01),
                  (6.783747177775549531e-02, -5.13898942407917780e-01, 8.551642850440059895e-01),
                  (-4.09811740153190107e-01, -6.51000036118079772e-01, -6.38947017057653221e-01),
                  (-4.35391566957816666e-01, 7.000228998401875069e-01, -5.66040743340399221e-01),
                  (-3.67275369838029608e-01, -9.06849357668608080e-01, -2.06719726214920074e-01),
                  (-7.24876157150795674e-01, 2.752762585295002173e-01, -6.31488351661607993e-01),
                  (7.999944917638109887e-01, -5.94809366836661745e-01, -7.88075521186358679e-02),
                  (-6.38660840611896762e-01, -7.30613237722315789e-01, -2.41488359002706959e-01),
                  (5.978837032381554284e-01, -9.09378167932159697e-02, 7.964078043810925989e-01),
                  (-9.29805720633694821e-01, -3.57570345696486047e-01, 8.720533100304592167e-02),
                  (4.682836169203208332e-01, -6.07627185281937154e-01, -6.41482390896564447e-01),
                  (7.309975652323822404e-01, 1.317784064962062851e-01, -6.69534921572915164e-01),
                  (-1.26449163913964823e-01, -9.90975538690972102e-01, -4.44757311530326326e-02),
                  (4.076648517742982869e-01, 5.753823407754260488e-01, 7.090448015123099745e-01),
                  (3.803260297965609382e-01, 9.217214425940288836e-01, 7.603744683751342825e-02),
                  (-7.25468075571367832e-01, 2.563257490022778362e-01, -6.38743439672923241e-01))).reshape(3,3,3,3)

def H_ext(T, Hmax, t, H0=0):
    while t > T:
       t -= T
    if t < T/4:
       H0 = (Hmax) *t/(T/4)
    elif t < 3*T/4:
       H0 = Hmax - 2 * Hmax * (t-(T/4))/(T/2)
    elif t < T:
       H0 = -Hmax + Hmax*(t-(3*T/4))/(T/4)
    return H0

Nx, Ny, Nz = 3, 3, 3
nx, ny, nz = N, N, N
dx, dy, dz = 900e-9 / (N-1), 900e-9 / (N-1), 900e-9 / (N-1)

mesh = Mesh(nx, ny, nz, dx, dy, dz)

m0 = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
m0[:, :, :, 2] = 1.

Ms = af.constant(0.0001/Constants.mu0, nx, ny, nz, dtype=af.Dtype.f64)
Ms[nx%Nx:,ny%Ny:,nz%Nz:] = 1.5/Constants.mu0
##Util.write_vti(Ms, dx, dy, dz, args.outdir + "Ms")

A = af.constant(0.0, nx, ny, nz, dtype=af.Dtype.f64)
A[nx%Nx:,ny%Ny:,nz%Nz:] = 10e-12
#Util.write_vti(A, dx, dy, dz, args.outdir + "A")

Ku = af.constant(0.0, nx, ny, nz, dtype=af.Dtype.f64)
Ku[nx%Nx:,ny%Ny:,nz%Nz:] = 8e3
#Util.write_vti(Ku, dx, dy, dz, args.outdir + "Ku")

#k_dir = np.loadtxt("k.txt").reshape(3,3,3,3)
Kdir = af.constant(1.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
Kdir[nx%Nx:,ny%Ny:,nz%Nz:,:] = af.from_ndarray(k_dir.repeat(nx//Nx,0).repeat(ny//Ny,1).repeat(nz//Nz,2))
#Util.write_vti(Kdir[:,:,:,0], dx, dy, dz, args.outdir + "Kdir_x")

# Create state object with timing
start = time.time()
state = State(mesh, Ms, m = m0)

if use_true_pbc:
    demag = DemagFieldPBC()
else:
    demag = DemagField(mesh, verbose = True)

exch = ExchangeFieldPBC(A,mesh)
aniso = UniaxialAnisotropyField(Ku, Kdir) 
zee = ExternalField(af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64))

llg = LLGIntegrator(alpha = 1, terms = [aniso, exch, zee, demag], rtol=1e-6, atol=1e-6)
with open(args.outdir + "m_gap%02.2fnm.dat" % (dx*1e9), "w") as fd:
    while state.t < 1e-9:
        llg.step(state)

    while state.t < 5/f/4 + 1e-9:
        t0 = state.t
        if use_true_pbc:
            H = H_ext(1/f, 0.1/Constants.mu0 , state.t-1e-9)
        else:
            H = H_ext(1/f, 2.5/Constants.mu0 , state.t-1e-9)
        zee.set_homogeneous_field(H, 0, 0)
        fd.write("%.15e, %.15e, %.15e, %.15e, %.15e, %.15e, %.15e, %.15e, \n" % (state.t, state.t-t0, H, 0, 0, state.mean_mx(), state.mean_my(), state.mean_mz()))
        fd.flush()
        llg.step(state)
