#!/usr/bin/python3
import arrayfire as af
import numpy as np
import magnumaf as maf
args = maf.parse()

def setup_regions(nx : int, ny : int, nz : int, dtype=af.Dtype.u32):
    """Define region array for multi-layer layout."""
    Regions = {"none" : 0, "skyrmion" : 1, "interlayer" : 2, "platinum" : 3}
    regions = af.constant(0, nx, ny, nz, dtype = dtype)
    regions[:, :, 0]  = Regions["skyrmion"]
    regions[:, :, 1]  = Regions["none"]
    regions[:, :, 2]  = Regions["none"]
    regions[:, :, 3]  = Regions["skyrmion"]
    regions[:, :, 4]  = Regions["none"]
    regions[:, :, 5]  = Regions["none"]
    regions[:, :, 6]  = Regions["skyrmion"]
    regions[:, :, 7]  = Regions["none"]
    regions[:, :, 8]  = Regions["none"]
    regions[:, :, 9]  = Regions["skyrmion"]
    regions[:, :, 10] = Regions["none"]
    regions[:, :, 11] = Regions["none"]
    regions[:, :, 12] = Regions["skyrmion"]
    regions[:, :, 13] = Regions["platinum"]
    regions[:, :, 14] = Regions["interlayer"]
    regions[:, :, 15] = Regions["interlayer"]
    regions[:, :, 16] = Regions["interlayer"]
    regions[:, :, 17] = Regions["interlayer"]
    regions[:, :, 18] = Regions["none"]
    regions[:, :, 19] = Regions["none"]
    regions[:, :, 20] = Regions["skyrmion"]
    regions[:, :, 21] = Regions["none"]
    regions[:, :, 22] = Regions["none"]
    regions[:, :, 23] = Regions["skyrmion"]
    regions[:, :, 24] = Regions["none"]
    regions[:, :, 25] = Regions["none"]
    regions[:, :, 26] = Regions["skyrmion"]
    regions[:, :, 27] = Regions["none"]
    regions[:, :, 28] = Regions["none"]
    regions[:, :, 29] = Regions["skyrmion"]
    regions[:, :, 30] = Regions["none"]
    regions[:, :, 31] = Regions["none"]
    regions[:, :, 32] = Regions["skyrmion"]
    return regions

def setup_init_skyrmion(dims_3d : [int], R = 0.25):
    """ set m [0,0,-1] for r <= R, else [0,0,1]"""
    nx, ny, nz = dims_3d[0], dims_3d[1], dims_3d[2]
    m = np.zeros((nx, ny, nz, 3))
    for ix in range (0, nx):
        for iy in range(0, ny):
            a= nx/2
            b= ny/2
            rx=ix-nx/2.
            ry=iy-ny/2.
            r = pow(rx, 2)/pow(a, 2)+pow(ry, 2)/pow(b, 2)
            if(r <= R):
                m[ix, iy, :, 2] = -1.
            else:
                m[ix, iy, :, 2] = 1.
    return af.from_ndarray(m)
    
def setup_m(regions):
    """Setup m, applying magnetization on geometry."""
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
