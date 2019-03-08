# based on: https://github.com/arrayfire/arrayfire-docker 
# image: magnum.af
# build: nvidia-docker build -t magnum.af -f Dockerfile --build-arg user="$UID" .
# run  : nvidia-docker run --rm -ti magnum.af /bin/bash
# test : nvidia-docker run --rm -t magnum.af ./tests/runall.sh .

FROM nvidia/cuda:9.2-devel-ubuntu18.04
MAINTAINER none

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
    make -j && \
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
             -DBUILD_CUDA=ON \
             -DBUILD_OPENCL=ON \
             -DBUILD_UNIFIED=ON \
             -DBUILD_GRAPHICS=${COMPILE_GRAPHICS} \
             -DBUILD_NONFREE=OFF \
             -DBUILD_EXAMPLES=ON \
             -DBUILD_TEST=ON \
             -DBUILD_DOCS=OFF \
             -DINSTALL_FORGE_DEV=ON \
             -DUSE_FREEIMAGE_STATIC=OFF && \
             # -DCOMPUTES_DETECTED_LIST="30;35;37;50;52;60" \
    make -j && make install && \
    mkdir -p ${AF_PATH} && ln -s /opt/arrayfire-3/* ${AF_PATH}/ && \
    echo "${AF_PATH}/lib" >> /etc/ld.so.conf.d/arrayfire.conf && \
    echo "/usr/local/cuda/nvvm/lib64" >> /etc/ld.so.conf.d/arrayfire.conf && \
    ldconfig && \
    cd ../.. && rm -rf arrayfire # saving >700MB

ENV ArrayFire_DIR=$ArrayFire_DIR:/opt/arrayfire/share/ArrayFire/cmake/ DEBIAN_FRONTEND=noninteractive DEBCONF_NONINTERACTIVE_SEEN=true

RUN apt update && \
    apt install -y python3-pip \
        libvtk7-dev && \
    pip3 install arrayfire cython numpy

# Instal Google Tests
RUN apt-get update && \
    apt-get install -y libgtest-dev && \
    cd /usr/src/gtest && \
    cmake CMakeLists.txt && \
    make && \
    cp *.a /usr/lib

# Setting user from build-arg with 999 as default
ARG user=999
RUN groupadd -g $user magnum.af.user && \
    useradd -r -u $user -g magnum.af.user magnum.af.user && \
    mkdir /home/magnum.af.user && \
    chown -R magnum.af.user /home/magnum.af.user
USER magnum.af.user

# Add magnum.af repository
COPY --chown=magnum.af.user . /home/magnum.af

WORKDIR /home/magnum.af

RUN scripts/magnum.af -vf -o build/main_empty scripts/main_empty.cpp && \
    ./tests/unit/cpp/maketests.sh . && \
    ./tests/integration/cpp/maketests.sh .
