#!/bin/bash
loops=1000
nnz=30
../magnum.af -f -p plot_timing.gpi timing.cpp "$HOME"/data_magnum.af/nonequi_demag/timing "$loops" "$nnz"
