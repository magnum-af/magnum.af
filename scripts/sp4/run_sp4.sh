#!/bin/bash

# calling this script's directory (to make it executable form everywhere)
cd $( dirname "${BASH_SOURCE[0]}" )

# running magnum.af
../magnum.af -f -p plot_sp4.gpi sp4.cpp $HOME/data_magnum.af/sp4/
