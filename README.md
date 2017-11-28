PTH-MAG: A general purpose magnetic simulation software
=====

## Main Features:
* Atomistic LLG Solver
* Micromagnetic LLG Solver
* String Method


## Prerequisites:
* A C++14 compiler, like gcc or clang
* [CMake](http://www.cmake.org) 3.0.0 or newer
* ArrayFire 3.0.1 or newer via. [pre-built binaries](http://arrayfire.com/download) or
  [source](https://github.com/arrayfire/arrayfire)

# Installation Guide:


## OpenCL Devices (e.g. AMD Graphics Cards):
* Installation of hardware-specific drivers:

http://support.amd.com/en-us/kb-articles/Pages/AMDGPU-PRO-Install.aspx
* AMD APP SDK 

http://developer.amd.com/amd-accelerated-parallel-processing-app-sdk/

* Fix linking error

Currently, AMDAPPSDK-3.0 installs with a broken symbolic link on 64-bit installs
which can be fixed by

going to directory
`/AMDAPPSDKROOT/lib/x86_64/` (usually `/opt/AMDAPPSDK-3.0/lib/x86_64`) 

and executing:

`$ sudo ln -sf sdk/libOpenCL.so.1 libOpenCL.so`

## Arrayfire 
* from binaries 

http://arrayfire.org/docs/installing.htm

* from source

https://github.com/arrayfire/arrayfire/wiki/Build-Instructions-for-Linux

## Arrayfire-Python


* Install arrayfire-python bindings by

`$ pip install arrayfire`

Fore more details see https://github.com/arrayfire/arrayfire-python

## VTK:
### From Package Manager
* Check current versions:

`$ apt-cache search libvtk`

* Choose latest version with headers  (i.e. latest -dev version), e.g. 

`$ sudo apt install libvtk6-dev`

Currently, this comes with a broken link in Ubuntu 16.04. Fix by installing

`$ sudo apt-get install libproj-dev`

An other quick and dirty fix is

`$ sudo ln -s  /usr/lib/x86_64-linux-gnu/libproj.so.<your-version> /usr/lib/x86_64-linux-gnu/libproj.so`

with <your-version> being eg 9 (https://stackoverflow.com/questions/37369369/compiling-pcl-1-7-on-ubuntu-16-04-errors-in-cmake-generated-makefile)

Note: maybe python-vtk is also neccessary
`$ sudo apt install python-vtk`
### From source (not recommended)
* follow

https://www.vtk.org/Wiki/VTK/Configure_and_Build

* generate ccmake files with

`$ ccmake /...path_to_VTK/`
* install with

`$ sudo make install -jXX`  (where XX is number of desired threads)

Note: if you encounter the error  "CMake Error at Rendering/OpenGL/CMakeLists.txt:304 (message):
   X11_Xt_LIB could not be found.  Required for VTK X lib."
try

`$ sudo apt-get install libxt-dev`

(https://stackoverflow.com/questions/23528248/how-to-install-x11-xt-lib-when-configure-vtk)

## PTH-MAG:
Note: In new projects, set VTK_DIR by

`$ export VTK_DIR=/home/.../VTK-build:$VTK_DIR`

* in pth_mag folder run

`$ mkdir build && cd build`

`$ cmake -DVTK_DIR:PATH=/home/path_to_VTK_build/VKT-build  ..`

`$ make -jXX` (where XX is number of desired threads)

Note: In new projects, set VTK_DIR by 
`$ export VTK_DIR=/home/.../VTK-build:$VTK_DIR`
