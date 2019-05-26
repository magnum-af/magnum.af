magnum.af: A finite differences GPU-accelerated micromagnetic and atomistic simulation software
=====
# Physical Methods
* Micromagnetic Model
* Atomistic Spin Model
* Solvers for the time-dependent Landau–Lifshitz–Gilbert equation
* Micromagnetic and Atomistic Energy Minimization
* Stochastic Langevin Dynamics
* String Method for Energy Barrier Calculations


# Main Features
* Support for CUDA, OpenCL and CPU backends.
  * This enables a high degree of flexibility in terms of user hardware as the
    code runns on both Nvidia(R) and AMD(R) devices as well as on any x86 CPU
* C++ Project
  * Optimized for performance
* Python bindings 
  * For an easy user-interface

# Installation Guide
## Docker:
For GPU support build the image provided in the Dockerfile by running the following command in the project's root directory:

`$ nvidia-docker build -t magnum.af -f Dockerfile --build-arg user="$UID" .`

For CPU support only use:

`$ docker build -t magnum.af.cpu -f Dockerfile.cpu --build-arg user="$UID" .`

For running simulations, use the provided script in 'scripts/magnum.af.docker', e.g.:

`$ magnum.af.docker sp4.py`

## Build scipt:
Execute the provided installation script:

`$./scripts/install_magnum.af_environment.sh`

## Manual installation:

### Prerequisites:
* A C++11 compiler like gcc or clang
* [CMake](http://www.cmake.org) 3.0.0 or newer
* pip3

### NVIDIA driver and CUDA:
following linuxconfig.org [linuxconfig.org](https://linuxconfig.org/how-to-install-the-nvidia-drivers-on-ubuntu-18-04-bionic-beaver-linux)

`$ sudo ubuntu-drivers autoinstall`

add your user to the video group

`$ sudo usermod -a -G video $LOGNAME`

Note: 
if the driver version provided by the repo is not sufficient, use a ppa instead:

`$ sudo add-apt-repository ppa:graphics-drivers && sudo apt-get update`

And select the proper driver version, e.g.:

`$ sudo apt install nvidia-driver-418 nvidia-settings`


install CUDA with

`$ sudo apt install nvidia-cuda-toolkit`

### Install Arrayfire
 For version 3.6.2:
 
`$ wget https://arrayfire.s3.amazonaws.com/3.6.2/ArrayFire-v3.6.2_Linux_x86_64.sh .`

`$ chmod +x ArrayFire-v3.6.2_Linux_x86_64.sh`

`$ sudo ./ArrayFire-v3.6.2_Linux_x86_64.sh  --include-subdir --prefix=/opt`

`$ sudo pip3 install arrayfire`

### Install VTK
`$ sudo apt install libvtk7-dev`

### Install Cython
`$ sudo pip3 install cython`

--------------------

### OpenCL Devices (e.g. AMD Graphics Cards):
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

### Arrayfire 
* from binaries 

http://arrayfire.org/docs/installing.htm

* from source

https://github.com/arrayfire/arrayfire/wiki/Build-Instructions-for-Linux

### Arrayfire-Python


* Install arrayfire-python bindings by

`$ pip install arrayfire`

Fore more details see https://github.com/arrayfire/arrayfire-python

### VTK:
#### From Package Manager
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

Add the VTK_DIR environment variable in your .bashrc file and source it

add `export VTK_DIR=/path_to_VTK:$VTK_DIR`

`source .bashrc`

#### From source (not recommended)
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

###  Arrayfire Python and Cython
$ pip install arrayfire
$ pip install Cython

### magnum.af:
Note: In new projects, set VTK_DIR by

`$ export VTK_DIR=/home/.../VTK-build:$VTK_DIR`

* in magnum.af folder run

`$ mkdir build && cd build`

`$ cmake -DVTK_DIR:PATH=/home/path_to_VTK_build/VKT-build  ..`

`$ make -jXX` (where XX is number of desired threads)

Note: In new projects, set VTK_DIR by 
`$ export VTK_DIR=/home/.../VTK-build:$VTK_DIR`
