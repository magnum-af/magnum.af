#!/bin/bash
# script installing all needed packages to run magnum.af
# NOTE: experimantal, use with care and skip not needed parts

# essentials
sudo apt -y install cmake gcc ipython

## Prerequisite for magnum.af, i.e. mainly arrayfire and vtk
sudo apt -y install build-essential libfreeimage3 libfontconfig1 libglu1-mesa python-pip libvtk7-dev
# Optional: paraview
sudo apt install paraview

# Cython
sudo pip install Cython

# Arrayfire 3.6.1 (update version if available)
wget -P ~/Downloads https://arrayfire.s3.amazonaws.com/3.6.2/ArrayFire-v3.6.2_Linux_x86_64.sh
chmod +x ~/Downloads/ArrayFire-v3.6.2_Linux_x86_64.sh
sudo ~/Downloads/ArrayFire-v3.6.2_Linux_x86_64.sh --include-subdir --prefix=/opt
echo "export LD_LIBRARY_PATH=/opt/arrayfire/lib64:\$LD_LIBRARY_PATH" >> ~/.bashrc
sudo pip install arrayfire

# Install Gtest
sudo apt-get install cmake # install cmake
cd /usr/src/gtest
sudo cmake CMakeLists.txt
sudo make
sudo cp *.a /usr/lib # copy or symlink libgtest.a and libgtest_main.a to your /usr/lib folder

# Install AMD GPU DRIVER (For AMD GPUs only):
# from  https://einsteinathome.org/content/quick-guide-how-install-opencl-amd-gpus-linux-kubuntu-1804-and-similar-distro
# compare https://www.amd.com/en/support/kb/release-notes/amdgpu-installation
wget -P ~/Downloads https://drivers.amd.com/drivers/linux/amdgpu-pro-18.40-676022-ubuntu-18.04.tar.xz # If wget is not sucessfull:
# install link for GTX580 Driver: https://www.amd.com/en/support/graphics/radeon-500-series/radeon-rx-500-series/radeon-rx-580
cd ~/Downloads
tar -xJpf amdgpu-pro-18.40-676022-ubuntu-18.04.tar.xz
cd amdgpu-pro-18.40-676022-ubuntu-18.04/
./amdgpu-install -y
sudo reboot
# TODO # run this after reboot:
cd ~/Downloads/amdgpu-pro-18.40-676022-ubuntu-18.04/
sudo ./amdgpu-pro-install -y --opencl=legacy
sudo reboot
# Uninstall with:
# sudo amdgpu-uninstall
