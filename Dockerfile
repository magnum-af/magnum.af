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
