#!/bin/bash
loops=1000
nnz=20
../../magnum.af -p plot_timing.gpi -b opencl timing.cpp "$HOME"/data_magnum.af/nonequi_demag/timing/opencl "$loops" "$nnz"
