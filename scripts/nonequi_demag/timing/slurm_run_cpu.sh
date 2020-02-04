#!/bin/bash
#SBATCH --job-name=nonequi_timing
#SBATCH --ntasks=1
#SBATCH --mem=1024
#SBATCH --time=70:00
loops=1000
nnz=20
../../magnum.af -p plot_timing.gpi -b cpu timing.cpp "$HOME"/data_magnum.af/nonequi_demag/timing/nx256_ny256/slurm_v0/cpu "$loops" "$nnz"
