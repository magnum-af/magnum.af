# build: $ docker build -t magnum.af -f Dockerfile --build-arg user="$UID" .
# push:  $ docker push magnum.af
# test:  $ docker run --rm -t magnum.af ./scripts/runalltests.sh .
# run:   $ docker run --rm -ti magnum.af /bin/bash

FROM nvidia/cuda:10.1-devel-ubuntu18.04

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    git \
    cmake \
    vim \
    # af:
    # avoiding boost 'fatal error: pyconfig.h':
    python-dev \
    #libboost-all-dev \ # manual below
    libglfw3-dev \
    libglu1-mesa-dev \
    clinfo \
    libfreeimage-dev \
    libfftw3-dev \
    libfontconfig1-dev \
    libfreeimage-dev \
    liblapacke-dev \
    libopenblas-dev \
    ocl-icd-opencl-dev \
    opencl-headers \
    wget \
    # magaf:
    software-properties-common \
    libvtk7-dev \
    ipython3 \
    python3-pip \
    # magaf docu:
    doxygen \
    python-pydot \
    python-pydot-ng \
    graphviz \
    && \
    # Pumping g++ version to 9
    add-apt-repository ppa:ubuntu-toolchain-r/test && \
    apt update && \
    DEBIAN_FRONTEND=noninteractive apt install -y \
    g++-9 && \
    update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-9 9 && \
    update-alternatives --auto g++ \
    #g++ --version
    && \
    rm -rf /var/lib/apt/lists/*

# installing boost 1.73
# TODO move before apt installs if fails
RUN wget https://dl.bintray.com/boostorg/release/1.73.0/source/boost_1_73_0.tar.bz2 && \
    tar --bzip2 -xf boost_*.tar.bz2 && \
    cd boost_* && \
    ./bootstrap.sh --prefix=/usr/include/boost/ && \
    ./b2 && \
    ./b2 install && \
    cd ../ && rm -r boost_*

# building arrayfire
RUN git clone --recursive https://github.com/arrayfire/arrayfire.git -b v3.7.2 && \
    cd arrayfire && mkdir build && cd build && \
    cmake .. -DBOOST_ROOT=/usr/include/boost/ && \
    make -j && \
    make install && \
    echo usr/local/lib/ > /etc/ld.so.conf.d/arrayfire.conf && \
    ldconfig && \
    rm -rf ../../arrayfire

# Setting up symlinks for libcuda and OpenCL ICD
RUN ln -s /usr/local/cuda/lib64/stubs/libcuda.so /usr/lib/libcuda.so.1 && \
    ln -s /usr/lib/libcuda.so.1 /usr/lib/libcuda.so && \
    mkdir -p /etc/OpenCL/vendors && \
    echo "libnvidia-opencl.so.1" > /etc/OpenCL/vendors/nvidia.icd && \
    echo "/usr/local/nvidia/lib" >> /etc/ld.so.conf.d/nvidia.conf && \
    echo "/usr/local/nvidia/lib64" >> /etc/ld.so.conf.d/nvidia.conf
ENV PATH=/usr/local/nvidia/bin:/usr/local/cuda/bin:${PATH}

# Install pip packages and Google Test
RUN apt-get update && \
    apt-get install -y libgtest-dev \
    python3-setuptools && \
    pip3 install arrayfire cython numpy && \
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

# Add magnum.af repository
COPY --chown=magnum.af.user . /home/magnum.af/

# building magnum.af and docu
WORKDIR /home/magnum.af/
RUN (mkdir build && cd build && cmake .. && make -j && make install) && \
    doxygen .doxygen-config && \
    chmod -R 777 /home/magnum.af/

# set non-root user
USER magnum.af.user
ENV HOME=/home/magnum.af.user \
    PYTHONPATH=/home/magnum.af/build/python/ \
    LD_LIBRARY_PATH=usr/local/lib/
