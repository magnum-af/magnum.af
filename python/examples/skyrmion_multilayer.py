#!/usr/bin/python3
import arrayfire as af
import numpy as np
import magnumaf as maf
args = maf.parse()

# Define layers as regions
def setup_layers(nx : int):
    Regions = {"none" : 0, "skyrmion" : 1, "interlayer" : 2, "platinum" : 3}
    layers = [0] * 33
    layers [0] = Regions["skyrmion"]
    layers [1] = Regions["none"]
    layers [2] = Regions["none"]
    layers [3] = Regions["skyrmion"]
    layers [4] = Regions["none"]
    layers [5] = Regions["none"]
    layers [6] = Regions["skyrmion"]
    layers [7] = Regions["none"]
    layers [8] = Regions["none"]
    layers [9] = Regions["skyrmion"]
    layers[10] = Regions["none"]
    layers[11] = Regions["none"]
    layers[12] = Regions["skyrmion"]
    layers[13] = Regions["platinum"]
    layers[14] = Regions["interlayer"]
    layers[15] = Regions["interlayer"]
    layers[16] = Regions["interlayer"]
    layers[17] = Regions["interlayer"]
    layers[18] = Regions["none"]
    layers[19] = Regions["none"]
    layers[20] = Regions["skyrmion"]
    layers[21] = Regions["none"]
    layers[22] = Regions["none"]
    layers[23] = Regions["skyrmion"]
    layers[24] = Regions["none"]
    layers[25] = Regions["none"]
    layers[26] = Regions["skyrmion"]
    layers[27] = Regions["none"]
    layers[28] = Regions["none"]
    layers[29] = Regions["skyrmion"]
    layers[30] = Regions["none"]
    layers[31] = Regions["none"]
    layers[32] = Regions["skyrmion"]
    return layers

def setup_regions(nx : int, ny : int, nz : int, dtype=af.Dtype.u32):
    layers = setup_layers(nz)
    regions = af.constant(0, 1, 1, len(layers), dtype = dtype)
    for i in range(len(layers)):
        regions[:, :, i] = layers[i]
    return af.tile(regions, nx, ny)

def setup_init_skyrmion(dims_3d : [int], r = 0.25):
    nx, ny, nz = dims_3d[0], dims_3d[1], dims_3d[2]
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
    
def setup_m(regions):
    binary_region_as_vec = maf.make_nonzeros_ones(regions)
    return af.tile(binary_region_as_vec, 1, 1, 1, 3) * setup_init_skyrmion(regions.dims())

### Config ###
nx, ny, nz = 100, 100, 33
dx = 400e-9/nx;
dy = 400e-9/ny;
dz = 1e-9;
print(nx, ny, nz, dx, dy, dz)

Hz_in_T = 130e-3
Hz_in_Apm = Hz_in_T / maf.Constants.mu0

# Material Parameters
RKKY_value = 0.8e-3 * dz
RKKY_layers = [13, 14]

# lookup values, position encodes region index
# region:      0,       1,       2,        3
Ms_values = [0.0,  1371e3, 488.2e3, 488.2e3/100.]
A_values  = [0.0,  15e-12,   4e-12, 2 * 15e-12  ]
Ku_values = [0.0, 1.411e6, 486.6e3, 0.0         ]

# In case1 DMI=0.0 in interlayer, in cases 2 and 3, DMI=0.8e-3
case2_or_3 = True
if (case2_or_3):
    D_values = [0.0, -2.5e-3, 0.8e-3, 0.0]
else:
    D_values = [0.0, -2.5e-3, 0.0, 0.0]
### End Config ###

regions = setup_regions(nx, ny, nz);
Ms = maf.lookup(Ms_values, regions)
A  = maf.lookup(A_values,  regions)
Ku = maf.lookup(Ku_values, regions)
D  = maf.lookup(D_values,  regions)

RKKY = af.constant(0.0, nx, ny, nz, dtype=af.Dtype.f64)
for layer in RKKY_layers:
    RKKY[:, :, layer] = RKKY_value

mesh = maf.Mesh(nx, ny, nz, dx, dy, dz)
state = maf.State(mesh, Ms = Ms, m = setup_m(regions))
state.write_vti(args.outdir + "m0")

# Interactions
exc = maf.RKKYExchangeField(RKKY, A, mesh)
dmg = maf.DemagField(mesh, verbose = True, caching = True, nthreads = 8)
uni = maf.UniaxialAnisotropyField(Ku, [0, 0, 1])
dmi = maf.DmiField(D, [0, 0, -1])
ext_field = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
ext_field[:, :, :, 2] = Hz_in_Apm
ext = maf.ExternalField(ext_field)
llg = maf.LLGIntegrator(alpha = 1, terms = [exc, dmg, uni, dmi, ext])

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

state.write_vti(args.outdir + "m_relaxed_" + str(relax_time_in_s))
