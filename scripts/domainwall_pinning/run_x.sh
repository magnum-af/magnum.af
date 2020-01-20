#!/bin/bash
cd "$(dirname "${BASH_SOURCE[0]}")"
#../magnum.af -fs -p plot_x.sh domainwall_pinning_x.py "$HOME"/data/domanwall_pinning/x/

../magnum.af -fS -p plot_x.sh domainwall_pinning_x_conv.py "$HOME"/data_magnum.af/domanwall_pinning/x_conv/
../magnum.af -fS -p plot_x.sh domainwall_pinning_x_sparse.py "$HOME"/data_magnum.af/domanwall_pinning/x_sparse/

../magnum.af -S -p plot_x.sh domainwall_pinning_x_conv.py "$HOME"/data_magnum.af/domanwall_pinning/x_conv/t100ns
../magnum.af -S -p plot_x.sh domainwall_pinning_x_sparse.py "$HOME"/data_magnum.af/domanwall_pinning/x_sparse/t100ns
