#!/bin/bash
loops=1000
nnz=20
../../magnum.af -p plot_timing.gpi -b cpu timing.cpp "$HOME"/data_magnum.af/nonequi_demag/timing/cpu "$loops" "$nnz"
