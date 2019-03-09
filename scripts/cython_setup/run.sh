#!/bin/bash
cd "$(dirname "${BASH_SOURCE[0]}")"
./cleanup.sh
cp ../../src/magnum_af.pyx .
cp ../../src/magnum_af_decl.pxd .
python3 setup.py build_ext --inplace 
python3 test.py
