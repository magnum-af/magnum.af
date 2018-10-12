FROM phusion/baseimage:0.10.1
MAINTAINER PTH <paul.thomas.heistracher@univie.ac.at> 
#image: magnum.af.opencl
#build:  docker build -t magnum.af.opencl -f Dockerfile.opencl .
#run image: docker run -ti magnum.af.opencl /bin/bash
#run tests: docker run --device=/dev/nvidia3 --device=/dev/nvidiactl --device=/dev/nvidia-uvm --device=/dev/nvidia-uvm-tools -t magnum.af.opencl /magnum.af/tests/runall.sh /magnum.af

# A docker container with the Nvidia kernel module and CUDA drivers installed 

RUN apt-get update && apt-get install -q -y wget build-essential kmod

RUN cd /opt && \ 
    wget https://developer.nvidia.com/compute/cuda/9.0/Prod/local_installers/cuda_9.0.176_384.81_linux-run && \
    chmod +x cuda_9.0.176_384.81_linux-run  

RUN cd /opt && \
    mkdir nvidia_installers && \
    ./cuda_9.0.176_384.81_linux-run -extract=`pwd`/nvidia_installers

RUN cd /opt/nvidia_installers && \ 
    ./NVIDIA-Linux-x86_64-384.81.run -s --no-kernel-module 
    
RUN cd /opt/nvidia_installers && \ 
    ./cuda-linux.9.0.176-22781540.run -noprompt 

RUN cd /opt/nvidia_installers && \
    ./cuda-samples.9.0.176-22781540-linux.run -noprompt -cudaprefix=/usr/local/cuda-9.0/

# Ensure the CUDA libs and binaries are in the correct environment variables 
ENV LIBRARY_PATH=$LD_LIBRARY_PATH:opt/arrayfire/lib:/opt/arrayfire:/usr/local/cuda-9.0/lib64 
ENV PATH=$PATH:/usr/local/cuda-9.0/bin

## Installing ArrayFire


RUN apt-get update && apt-get install -y software-properties-common && \
    apt-get install -y --no-install-recommends \
        build-essential \
        clinfo \
        cmake \
        git \
        libboost-all-dev \
        libfftw3-dev \
        libfontconfig1-dev \
        libfreeimage-dev \
        liblapack-dev \
        liblapacke-dev \
        libopenblas-dev \
        ocl-icd-opencl-dev \
        opencl-headers \
        wget \
        xorg-dev && \
    rm -rf /var/lib/apt/lists/*

# Setting up symlinks for libcuda and OpenCL ICD
RUN ln -s /usr/local/cuda/lib64/stubs/libcuda.so /usr/lib/libcuda.so.1 && \
    ln -s /usr/lib/libcuda.so.1 /usr/lib/libcuda.so && \
    mkdir -p /etc/OpenCL/vendors && \
    echo "libnvidia-opencl.so.1" > /etc/OpenCL/vendors/nvidia.icd && \
    echo "/usr/local/nvidia/lib" >> /etc/ld.so.conf.d/nvidia.conf && \
    echo "/usr/local/nvidia/lib64" >> /etc/ld.so.conf.d/nvidia.conf
ENV PATH=/usr/local/nvidia/bin:/usr/local/cuda/bin:${PATH}

WORKDIR /root

# Build GLFW from source
RUN git clone https://github.com/glfw/glfw.git && \
    cd glfw && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=/usr .. && \
    make -j16 && \
    make install


# AF_DISABLE_GRAPHICS - Environment variable to disable graphics at
# runtime due to lack of graphics support by docker - visit
# http://arrayfire.org/docs/configuring_environment.htm#af_disable_graphics
# for more information
ENV AF_PATH=/opt/arrayfire AF_DISABLE_GRAPHICS=1
ARG COMPILE_GRAPHICS=OFF
RUN git clone --recursive https://github.com/arrayfire/arrayfire.git -b master && \
    cd arrayfire && mkdir build && cd build && \
    cmake .. -DCMAKE_INSTALL_PREFIX=/opt/arrayfire-3 \
             -DCMAKE_BUILD_TYPE=Release \
             -DBUILD_CPU=ON \
             -DBUILD_CUDA=OFF \
             -DBUILD_OPENCL=ON \
             -DBUILD_UNIFIED=OFF \
             -DBUILD_GRAPHICS=${COMPILE_GRAPHICS} \
             -DBUILD_NONFREE=OFF \
             -DBUILD_EXAMPLES=OFF \
             -DBUILD_TEST=OFF \
             -DBUILD_DOCS=OFF \
             -DINSTALL_FORGE_DEV=ON \
             -DUSE_FREEIMAGE_STATIC=OFF && \
             # -DCOMPUTES_DETECTED_LIST="30;35;37;50;52;60" \
    make -j16 && make install && \
    mkdir -p ${AF_PATH} && ln -s /opt/arrayfire-3/* ${AF_PATH}/ && \
    echo "${AF_PATH}/lib" >> /etc/ld.so.conf.d/arrayfire.conf && \
    echo "/usr/local/cuda/nvvm/lib64" >> /etc/ld.so.conf.d/arrayfire.conf && \
    ldconfig

WORKDIR /root/arrayfire

# Installing VTK
RUN apt-get update && apt-get -y install libvtk6-dev
RUN apt-get -y install libproj-dev
ENV VTK_DIR=$VTK_DIR:/usr/lib/tcltk/vtk-6.2

RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
RUN python get-pip.py
RUN pip install Cython
RUN pip install arrayfire
RUN pip install numpy
# Install GFLW
RUN apt-get -y install libglfw3-dev

# Instal Google Tests
RUN apt-get install libgtest-dev && \
    cd /usr/src/gtest && \
    cmake CMakeLists.txt && \
    make && \
    cp *.a /usr/lib

# Add magnum.af repository
RUN mkdir /tmp/magnum.af

COPY . /tmp/magnum.af

#Setting variable for cmake, shout be set previously by "echo /opt/arrayfire/lib > /etc/ld.so.conf.d/arrayfire.conf && \ ldconfig" but is not working properly
ENV ArrayFire_DIR=$ArrayFire_DIR:/opt/arrayfire/share/ArrayFire/cmake/

RUN cd /tmp/magnum.af && \
    cp scripts/sp4/main_sp4.cpp src/ && \
    rm -rf build && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j12 && \
    cd .. && \
    ./tests/unit/cpp/maketests.sh . && \
    ./tests/integration/cpp/maketests.sh . 

RUN mkdir /magnum.af && mkdir /magnum.af/build && mkdir /magnum.af/build/src && \
    mv /tmp/magnum.af/build/src/magnum_af.so /magnum.af/build/src/ && \
    mv /tmp/magnum.af/tests /magnum.af && \
    rm -rf /tmp/magnum.af










## Installing ArrayFire
#RUN cd /opt && \
#    wget --no-check-certificate https://go.pardot.com/l/37882/2015-03-19/lw515 
#RUN cd /opt && \
#    chmod +x lw515
#
#RUN apt -y install build-essential cmake cmake-curses-gui
#
#RUN cd /opt && \
#    ./lw515 --include-subdir --prefix=/opt
#
## Adding Arrayfire libraries via ldconfig #TODO not working as settig ArrayFire_DIR is necessary
#RUN cd /opt && \
#    echo /opt/arrayfire/lib > /etc/ld.so.conf.d/arrayfire.conf && \
#    ldconfig
#
#RUN apt -y install build-essential libfreeimage3 libfontconfig1 libglu1-mesa
#
## Testing Arrayfire
#RUN mkdir todel 
#RUN cp -r /opt/arrayfire/share/ArrayFire/examples /todel/examples
#ENV ArrayFire_DIR=$ArrayFire_DIR:/opt/arrayfire/share/ArrayFire/cmake/
#RUN cd /todel/examples && \
#    ls && \
#    mkdir build && \
#    cd build
#
#RUN cd /todel/examples/build && \
#    cmake -DASSETS_DIR:PATH=/todel .. && \
#    make
#
#RUN cd /todel/examples/build && \
#    ./helloworld/helloworld_cpu  
#    #./helloworld/helloworld_opencl && \
#    #./helloworld/helloworld_cuda 
#
## Installing VTK
#RUN apt -y install libvtk6-dev
#RUN apt-get -y install libproj-dev
#ENV VTK_DIR=$VTK_DIR:/usr/lib/tcltk/vtk-6.2
#
#RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
#RUN python get-pip.py
#RUN pip install Cython
#RUN pip install arrayfire
## Install GFLW
#RUN apt-get -y install libglfw3-dev
#
## Add magnum.af repository
#RUN mkdir /tmp/magnum.af
#
#COPY . /tmp/magnum.af
#
#RUN cd /tmp/magnum.af && \
#    cp scripts/sp4/main_sp4.cpp src/ && \
#    mkdir build && \
#    cd build && \
#    cmake .. && \
#    make -j12

#TODO
#NOT WORKING: CUDA BACKEND ON GTO CONTAINER!!!
#AS a result, python gets stuck when loading af

#TODO
# Cleanup
#RUN rm -r /tmp/magnum.af

#NOTES:
#on GTO:
# docker build -t <docker_myname> .
# docker run -ti --device=/dev/nvidia0 --device=/dev/nvidiactl --device=/dev/nvidia-uvm --device=/dev/nvidia-uvm-tools --device=/dev/nvidia1 --device=/dev/nvidia2 --device=/dev/nvidia3  paul /bin/bash




########################################
#ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/arrayfire/lib 
#ENV ArrayFire_DIR=$ArrayFire_DIR:/opt/arrayfire/include 
#ENV ArrayFire_DIR=$ArrayFire_DIR:/opt/arrayfire/lib 
########################################
## Install VTK
#RUN cd /opt && \
#    git clone git://vtk.org/VTK.git VTK 
#    #cd VTK  && \
#    #git checkout tags/v8.1.0.rc3 -b v8.1.0.rc3
#
#RUN  cd /opt && \
#    mkdir VTK-build && \
#    cd VTK-build && \
#    cmake -DCMAKE_BUILD_TYPE:STRING=Release ../VTK && \
#    make -j12

