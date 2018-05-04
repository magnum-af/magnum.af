FROM phusion/baseimage:0.10.1
MAINTAINER PTH <paul.thomas.heistracher@univie.ac.at> 

# A docker container with the Nvidia kernel module and CUDA drivers installed 

RUN apt-get update && apt-get install -q -y wget build-essential kmod

RUN cd /opt && \ 
    wget https://developer.nvidia.com/compute/cuda/8.0/Prod2/local_installers/cuda_8.0.61_375.26_linux-run && \
    chmod +x cuda_8.0.61_375.26_linux-run  

RUN cd /opt && \
    mkdir nvidia_installers && \
    ./cuda_8.0.61_375.26_linux-run -extract=`pwd`/nvidia_installers

RUN cd /opt/nvidia_installers && \ 
    ./NVIDIA-Linux-x86_64-375.26.run -s --no-kernel-module 
    
RUN cd /opt/nvidia_installers && \ 
    ./cuda-linux64-rel-8.0.61-21551265.run -noprompt 
    
RUN cd /opt/nvidia_installers && \
    ./cuda-samples-linux-8.0.61-21551265.run -noprompt -cudaprefix=/usr/local/cuda-8.0/

# Ensure the CUDA libs and binaries are in the correct environment variables 
ENV LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-8.0/lib64 
ENV PATH=$PATH:/usr/local/cuda-8.0/bin

# Installing ArrayFire
RUN cd /opt && \
    wget --no-check-certificate https://go.pardot.com/l/37882/2015-03-19/lw515 
RUN cd /opt && \
    chmod +x lw515

RUN apt -y install build-essential cmake cmake-curses-gui

RUN cd /opt && \
    ./lw515 --include-subdir --prefix=/opt

# Adding Arrayfire libraries via ldconfig #TODO not working as settig ArrayFire_DIR is necessary
RUN cd /opt && \
    echo /opt/arrayfire/lib > /etc/ld.so.conf.d/arrayfire.conf && \
    ldconfig

RUN apt -y install build-essential libfreeimage3 libfontconfig1 libglu1-mesa

# Testing Arrayfire
RUN mkdir todel 
RUN cp -r /opt/arrayfire/share/ArrayFire/examples /todel/examples
ENV ArrayFire_DIR=$ArrayFire_DIR:/opt/arrayfire/share/ArrayFire/cmake/
RUN cd /todel/examples && \
    ls && \
    mkdir build && \
    cd build

RUN cd /todel/examples/build && \
    cmake -DASSETS_DIR:PATH=/todel .. && \
    make

RUN cd /todel/examples/build && \
    ./helloworld/helloworld_cpu  
    #./helloworld/helloworld_opencl && \
    #./helloworld/helloworld_cuda 

# Installing VTK
RUN apt -y install libvtk6-dev
RUN apt-get -y install libproj-dev
ENV VTK_DIR=$VTK_DIR:/usr/lib/tcltk/vtk-6.2

RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
RUN python get-pip.py
RUN pip install Cython
RUN pip install arrayfire
# Install GFLW
RUN apt-get -y install libglfw3-dev

# Add pth-mag repository
RUN mkdir /tmp/pth-mag

COPY . /tmp/pth-mag

RUN cd /tmp/pth-mag && \
    cp scripts/sp4/main_sp4.cpp src/ && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j12

#TODO
#NOT WORKING: CUDA BACKEND ON GTO CONTAINER!!!
#AS a result, python gets stuck when loading af

#TODO
# Cleanup
#RUN rm -r /tmp/pth-mag

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

