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

* fix linking error

"On 64-bit installs, AMDAPPSDK-3.0 installs with a broken symbolic link which 
can be fixed by going to $AMDAPPSDKROOT/lib/x86_64/ (usually $/opt/AMDAPPSDK-3.0/lib/x86_64) 

by  executing:

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
* follow

https://www.vtk.org/Wiki/VTK/Configure_and_Build

* generate ccmake files with

`$ ccmake /...path_to_VTK/`
* install with

`$ sudo make install -jXX`  (where XX is number of threads)

Note: if you encounter the error  "CMake Error at Rendering/OpenGL/CMakeLists.txt:304 (message):
   X11_Xt_LIB could not be found.  Required for VTK X lib."
try

`$ sudo apt-get install libxt-dev`

(https://stackoverflow.com/questions/23528248/how-to-install-x11-xt-lib-when-configure-vtk)

## PTH-MAG:
* Note in new projects, set VTK_DIR by

`$ export VTK_DIR=/home/.../VTK-build:$VTK_DIR`

* in pth_mag folder

`$ mkdir build && cd build`

`$ cmake -DVTK_DIR:PATH=/home/path_to_VTK_build/VKT-build  ..`


Note: In new projects, set VTK_DIR by 
`$ export VTK_DIR=/home/.../VTK-build:$VTK_DIR`


# TODEL
## Building this project

### Linux and OSX
To build this project simply:

```
git clone https://github.com/arrayfire/arrayfire-project-templates.git
cd arrayfire-project-templates/CMake/build
cmake ..
This is a base project which demonstrates how to properly configure CMake
to build a project that links with the [ArrayFire](http://arrayfire.com)
high performance parallel computing library. The CMake configuration files
in this directory will automatically find ArrayFire and build CPU, CUDA, and
OpenCL backends if support for all three backends is present on your machine.

make
```

NOTE: If you have installed ArrayFire to a non-standard location, you will
need to specify the full path to the `share/ArrayFire/cmake` directory. For
example, if you have installed ArrayFire to `/opt`, then the `cmake` command
above will be replaced with:

```
cmake -DArrayFire_DIR=/opt/share/ArrayFire/cmake ..
```

### Windows
Download the project and files into your working directory. For simplicity,
lets assume the project is at `C:\workspace\arrayfire-cmake`.

Open the CMake GUI.

In the source directory field, enter the working directory
`C:\workspace\arrayfire-cmake`. In the build directory field, enter the build
directory `C:\workspace\arrayfire-cmake\build`.

Click `Configure`.

If there are no errors, click `Generate` in the CMake GUI.

The ArrayFire installer creates files and registries that allow CMake to
detect it automatically. When configure is clicked, ArrayFire along with it's
dependencies are detected and the available backends are set to available.

When not using the installer, the ArrayFire directory will have to be set
using the `-DArrayFire_DIR` CMake variable.

Now, the solution files will be created in the build directory. There will be
different projects created for CPU, CUDA and OpenCL. Setting one of these
projects as start project and then running it will run the example with the
specified backend.
