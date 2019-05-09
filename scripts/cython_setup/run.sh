#!/bin/bash
pythonv=python3
cd "$(dirname "${BASH_SOURCE[0]}")"
./cleanup.sh
cp ../../src/magnumaf.pyx .
cp ../../src/magnumaf_decl.pxd .
"$pythonv" setup.py build_ext --inplace && "$pythonv" test.py
