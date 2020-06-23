#!/bin/bash
#../../magnum.af -p plot_inf_saf.gpi -f inf_saf.cpp $HOME/data_magnum.af/inf_saf/testrun/

# prev in examples dir:
#(cd ../build/ && make -j) && cd ../bin && ./inf_saf_cpu run_inf_saf/ && cat run_inf_saf/m.dat
#gnuplot plot_inf_saf.gpi

Hmax_mT="0.100"
for RKKY in $(LC_ALL=C seq -0 -0.1 -0.8); do
    echo $RKKY
    if ( -n "$HOME/data_magnum.af/inf_saf/run1/RKKYmT$RKKY/" ); then
    ../../magnum.af -p plot_inf_saf.gpi -f inf_saf.cpp $HOME/data_magnum.af/inf_saf/run1/RKKYmT$RKKY/ "$Hmax_mT" "$RKKY"
    fi
done
