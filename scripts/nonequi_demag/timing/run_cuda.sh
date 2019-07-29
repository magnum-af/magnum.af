#!/bin/bash
loops=1000
nnz=20
../../magnum.af -f -p plot_timing.gpi -b cuda timing.cpp "$HOME"/data_magnum.af/nonequi_demag/timing/cuda "$loops" "$nnz"
