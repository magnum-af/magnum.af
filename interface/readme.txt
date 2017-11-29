https://studywolf.wordpress.com/2012/09/14/cython-journey-part-1/

python setup.py build_ext -i

Troubleshooting:

* Forgetting to mention a file in sources in setup.py results in error like:

ImportError: /home/pth/git/pth-mag/interface/pth_mag.so: undefined symbol: _ZN4MeshC1Eiiiddd

Include respective file to solve

* Core dump at startup often results when choosing the wrong library in setup.py for desired backend (cpu, opencl, cuda)

E.g.: ibraries = ["afopencl"] 

Include correct library 

