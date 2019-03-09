#!/bin/bash
pythonv=python3
cd "$(dirname "${BASH_SOURCE[0]}")"
./cleanup.sh
cp ../../src/magnum_af.pyx .
cp ../../src/magnum_af_decl.pxd .
"$pythonv" setup.py build_ext --inplace && "$pythonv" test.py
