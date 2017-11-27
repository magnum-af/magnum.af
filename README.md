PTH-MAG: A general purpose magnetic simulation software
=====

## Main Feature:
* Atomistic LLG Solver
* Micromagnetic LLG Solver
* String Method

## Installation Guide: (PRELIMINARY)


* AMD Graphics Card:
Installation of hardware-specific drivers:
$ http://support.amd.com/en-us/kb-articles/Pages/AMDGPU-PRO-Install.aspx
* AMD APP SDK 
http://developer.amd.com/amd-accelerated-parallel-processing-app-sdk/

* VTK:
$ cmake -DVTK_DIR:PATH=/home/paul/Programs/VKT-build ..

Note: In new projects, set VTK_DIR by 
$ export VTK_DIR=/home/.../VTK-build:$VTK_DIR

## Prerequisites:
* A C++14 compiler, like gcc or clang
* [CMake](http://www.cmake.org) 3.0.0 or newer
* ArrayFire 3.0.1 or newer via. [pre-built binaries](http://arrayfire.com/download) or
  [source](https://github.com/arrayfire/arrayfire)

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
