magnum.af: A finite differences GPU-accelerated micromagnetic and atomistic simulation software
=====
# Physical Methods
* Micromagnetic Model
* Atomistic Spin Model
* Solvers for the time-dependent Landau–Lifshitz–Gilbert equation
* Micromagnetic and Atomistic Energy Minimization
* Stochastic Langevin Dynamics
* String Method for Energy Barrier Calculations


# Architecture
* Support for CUDA, OpenCL and CPU backends.
  * This enables a high degree of flexibility in terms of user hardware as the
    code runns on both Nvidia(R) and AMD(R) devices as well as on any x86 CPU
* C++ Project
  * Core functionality including all interactions and solvers
* Python bindings
  * Script-style user interface, compatible with numpy input arrays (see python/examples for more)

# Installation Guide
## Docker:
For GPU support build the image provided in the Dockerfile by running the following command in the project's root directory:

`$ nvidia-docker build -t magnum.af -f Dockerfile --build-arg user="$UID" .`

For CPU support only use:

`$ docker build -t magnum.af.cpu -f Dockerfile.cpu --build-arg user="$UID" .`

For running simulations, use the provided script in 'scripts/magnum.af.docker', e.g.:

`$ magnum.af.docker sp4.py`

## Local installation:

### GPU Driver Installation

Make sure your user is part of the video group using `$groups | grep video`. 

If not, add your user to the video group:

`$ sudo usermod -a -G video $LOGNAME`

perform a logout-login for the changes to have effect.

In the following, choose either NVIDIA or AMD depending on your hardware:

#### a) NVIDIA: nvidia-driver and CUDA:
following [linuxconfig.org](https://linuxconfig.org/how-to-install-the-nvidia-drivers-on-ubuntu-18-04-bionic-beaver-linux):

`$ sudo ubuntu-drivers autoinstall`

Note:
if the driver version provided by the repo is not sufficient, use a ppa instead:

`$ sudo add-apt-repository ppa:graphics-drivers && sudo apt-get update`

And select the proper driver version, e.g.:

`$ sudo apt install nvidia-driver-418 nvidia-settings`


install CUDA with

`$ sudo apt install nvidia-cuda-toolkit`

#### b) AMD: AMDGPU driver and OpenCL
Follow [linuxconfig.org](https://linuxconfig.org/how-to-install-the-latest-amd-radeon-drivers-on-ubuntu-18-04-bionic-beaver-linux)

On the tested system, the install script needed to be invoked with the option `--opencl=legacy`
(tested with driver: amdgpu-pro-19.20-812932-ubuntu-18.04 for Radeon 580 series GPU):

`$ ./amdgpu-pro-install --opencl=legacy -y`

for cmake to find OpenCl run:

`$ sudo apt install ocl-icd-opencl-dev`

### Install Arrayfire
 For version 3.6.2:

`$ wget https://arrayfire.s3.amazonaws.com/3.6.2/ArrayFire-v3.6.2_Linux_x86_64.sh .`

`$ chmod +x ArrayFire-v3.6.2_Linux_x86_64.sh`

`$ sudo ./ArrayFire-v3.6.2_Linux_x86_64.sh  --include-subdir --prefix=/opt`

For the loader to find the arrayfire lib use either:

`$ echo /opt/arrayfire/lib > /etc/ld.so.conf.d/arrayfire.conf`

`$ sudo ldconfig`

Or add the following line to your .bashrc file:

`$ echo 'export LD_LIBRARY_PATH=/opt/arrayfire/lib:$LD_LIBRARY_PATH' >> .bashrc`

### Install VTK
`$ sudo apt install libvtk7-dev`

### Install Python Packages
`$ sudo pip3 install cython`

`$ sudo pip3 install arrayfire`

### Build magnum.af
Prerequisites:
* A C++11 compiler like gcc or clang
* [CMake](http://www.cmake.org) 2.8 or newer

`$ mkdir build && cd build && cmake .. && make -j`

Install to system directories using:

`$ make install`

For python to find the shared library, add the install directory to the PYTHONPATH:

`$ echo 'export PYTHONPATH=/usr/local/lib/:$PYTHONPATH' >> .bashrc`

or alternatively, add a symbolic link to a folder in the PYTHONPATH, e.g. for python 3.6
$ `sudo ln -s /usr/local/lib/magnumaf.so /usr/lib/python3.6/
