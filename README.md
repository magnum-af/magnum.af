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
* [CMake](http://www.cmake.org) 2.8 or newer
* pip3

### GPU Driver Installation (choose either NVIDIA or AMD)
#### NVIDIA: nvidia-driver and CUDA:
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

#### AMD: AMDGPU driver and OpenCL
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

for cmake to find OpenCl run:

`$ sudo apt install ocl-icd-opencl-dev`

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
