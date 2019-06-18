#!/bin/bash

# calling this script's directory (to make it executable form everywhere)
cd $( dirname "${BASH_SOURCE[0]}" )

# running magnum.af
../../magnum.af -f -p plot_z_vs_x_sp4.gpi z_vs_x_sp4.cpp $PWD/output_z_vs_x_sp4
