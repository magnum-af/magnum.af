#!/usr/bin/env python3

from PYBINDTESTLIB import greet
print(greet())


from PYBINDTESTLIB import Mesh
mesh = Mesh(1, 2, 3, 0.1, 0.2, 0.3)
print(mesh.nx, mesh.ny, mesh.nz, mesh.dx, mesh.dy, mesh.dz, mesh.V)
