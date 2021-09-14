# build: $ docker build -t magnum.af -f Dockerfile --build-arg user="$UID" .
# test:  $ docker run --rm --gpus all -t magnum.af ./bash/runalltests.sh .
# run:   $ docker run --rm --gpus all -ti magnum.af /bin/bash

FROM nvidia/cuda:11.4.1-devel-ubuntu20.04
MAINTAINER Paul Heistracher <paul.thomas.heistracher@univie.ac.at>

# Important note for the following RUN command:
#   Arrayfire packages must be in seperate 'apt install' call!
#   Otherwise causing build warnings/errors:
#     '/usr/bin/ld: warning: libcuda.so.1, needed by /usr/lib/libafcuda.so.3.8.0, not found'
#     '/usr/bin/ld: /usr/lib/libafcuda.so.3.8.0: undefined reference to `cu<...>'
#   Somehow 'arrayfire-cuda3-dev' and 'arrayfire' cannot be in same 'apt install' as other packages.
#   Maybe some path becomes overwritten or not set otherwise?
# Note: Package 'ocl-icd-opencl-dev' needed for cmake to find OpenCL

RUN apt update && apt install -y gnupg curl && \
    curl -sSL "https://repo.arrayfire.com/GPG-PUB-KEY-ARRAYFIRE-2020.PUB" | apt-key add - && \
    echo "deb [arch=amd64] https://repo.arrayfire.com/ubuntu focal main" | tee /etc/apt/sources.list.d/arrayfire.list && \
    apt update && DEBIAN_FRONTEND=noninteractive apt install -y --no-install-recommends \
    build-essential \
    git \
    cmake \
    gnuplot \
    libboost-program-options-dev \
    libvtk7-dev \
    libgmock-dev \
    ocl-icd-opencl-dev \
    python3-pip && \
    pip3 install arrayfire numpy && \
    apt update && apt install -y \
    arrayfire-cuda3-dev \
    arrayfire

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
RUN chmod -R 777 /home/magnum.af/ && \
    (mkdir build && cd build && cmake .. && make -j && make install)

# set non-root user
USER magnum.af.user

ENV HOME=/home/magnum.af.user/ \
    PYTHONPATH=/usr/local/lib/
