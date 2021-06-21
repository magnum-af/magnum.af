#!/usr/bin/python3
import arrayfire as af
from magnumaf import *
from enum import Enum
args = parse()

Regions = {"vacuum" : 0, "skyrmion" : 1, "interlayer" : 2, "platinum" : 3}

def setup_regions(nx, ny, nz):
    regions = af.constant(0, 1, 1, nz, dtype=af.Dtype.u32)
    regions[:, :, 0]  = Regions["skyrmion"]
    regions[:, :, 3]  = Regions["skyrmion"]
    regions[:, :, 6]  = Regions["skyrmion"]
    regions[:, :, 9]  = Regions["skyrmion"]
    regions[:, :, 12] = Regions["skyrmion"]

    regions[:, :, 13] = Regions["platinum"]

    regions[:, :, 14] = Regions["interlayer"]
    regions[:, :, 15] = Regions["interlayer"]
    regions[:, :, 16] = Regions["interlayer"]
    regions[:, :, 17] = Regions["interlayer"]

    regions[:, :, 20] = Regions["skyrmion"]
    regions[:, :, 23] = Regions["skyrmion"]
    regions[:, :, 26] = Regions["skyrmion"]
    regions[:, :, 29] = Regions["skyrmion"]
    regions[:, :, 32] = Regions["skyrmion"]

    return af.tile(regions, nx, ny)

def setup_init_skyrmion(dims, r = 0.5):
    nx, ny, nz = dims[0], dims[1], dims[2]
    m = np.zeros((nx, ny, nz, 3))
    for ix in range (0, nx):
        for iy in range(0, ny):
            a= nx/2
            b= ny/2
            rx=ix-nx/2.
            ry=iy-ny/2.
            rr = pow(rx, 2)/pow(a, 2)+pow(ry, 2)/pow(b, 2)
            if(rr <= r):
                m[ix, iy, :, :] =-1.
            else:
                m[ix, iy, :, :] = 1.
    return af.from_ndarray(m)

    
def one_if_nonzero_else_zero(array):
    return af.iszero(af.iszero(array))

def setup_m(regions):
    binary_region_as_vec = one_if_nonzero_else_zero(regions)
    return af.tile(binary_region_as_vec, 1, 1, 1, 3) * setup_init_skyrmion(regions.dims(), 0.25)

def do_lookup(values, index, as_type = af.Dtype.f64):
    values_array = af.array.Array(values).as_type(as_type)
    look = af.lookup(values_array, af.flat(index))
    nx, ny, nz = index.dims()
    return af.moddims(look, nx, ny, nz)

### Config ###
nx, ny, nz = 100, 100, 33

regions = setup_regions(nx, ny, nz);
m0 = setup_m(regions)

dx = 400e-9/nx;
dy = 400e-9/ny;
dz = 1e-9;
print(nx, ny, nz, dx, dy, dz)

Hz_in_T = 130e-3
Hz_in_Apm = Hz_in_T / Constants.mu0

### Material Parameters ###
# position encode index, [0] is region 0, [1] is region 1, ...
Ms_values = [0.0, 1371e3, 488.2e3, 488.2e3/100.]
A_values  = [0.0, 15e-12, 4e-12, 2 * 15e-12]
Ku_values = [0.0, 1.411e6, 486.6e3, 0.0]

# In case1 DMI=0.0 in interlayer, in cases 2 and 3, DMI=0.8e-3
case2_or_3 = True
if (case2_or_3):
    D_values = [0.0, -2.5e-3, 0.8e-3, 0.0]
else:
    D_values = [0.0, -2.5e-3, 0.0, 0.0]

Ms = do_lookup(Ms_values, regions)
A = do_lookup(A_values, regions)
Ku = do_lookup(Ku_values, regions)
D = do_lookup(D_values, regions)

RKKY_value = 0.8e-3 * dz
RKKY = af.constant(0.0, nx, ny, nz, dtype=af.Dtype.f64)
RKKY[:, :, 13] = RKKY_value
RKKY[:, :, 14] = RKKY_value

# print("Ms", Ms)
# print("Ku", Ku)

mesh = Mesh(nx, ny, nz, dx, dy, dz)
state = State(mesh, Ms = Ms, m = m0)
state.write_vti(args.outdir + "m0")

# Interactions
exc = RKKYExchangeField(RKKY, A, mesh)
dmg = DemagField(mesh, verbose = True, caching = False, nthreads = 8)
uni = UniaxialAnisotropyField(Ku, [0, 0, 1])
dmi = DmiField(D, [0, 0, -1])
ext_field = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
ext_field[:, :, :, 2] = Hz_in_Apm
ext = ExternalField(ext_field)
llg = LLGIntegrator(alpha = 1, terms = [exc, dmg, uni, dmi, ext])

relax_time_in_s = 50e-9
stream = open(args.outdir + "m.dat", "w", buffering = 1)
vti_every_nth_step = 1000
while state.t < relax_time_in_s:
    llg.step(state)
    mx, my, mz = state.mean_m()
    Mx, My, Mz = state.mean_M()
    print(llg.accumulated_steps, state.t, mx, my, mz, Mx, My, Mz, flush = True)
    stream.write("%e, %e, %e, %e\n" %(state.t, mx, my, mz))
    if(llg.accumulated_steps % vti_every_nth_step == 0):
        state.write_vti(args.outdir + "m_step_" + str(llg.accumulated_steps))
