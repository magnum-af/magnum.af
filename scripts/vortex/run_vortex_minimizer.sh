#!/bin/bash
cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # cd this scripts dir
../magnum.af -f -p plot_hysteresis.sh main_vortex_minimizer.cpp $HOME/data_magnum.af/vortex/hysteresis_minimizer
