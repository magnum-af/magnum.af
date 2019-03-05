import sys
import os
import arrayfire as af
import numpy as np
from magnum_af import *
import time

print ("The arguments are: " , str(sys.argv))
filepath = sys.argv[1]

x = 1000e-9
y = 1000e-9
z = 80e-9

nx = 250
ny = 250
nz = 1

disk, n_cells = Util.disk(nx, ny, nz, 0)
print (disk)
