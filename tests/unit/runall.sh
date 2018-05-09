#!/bin/bash

#usage: ./runall /path/to/gitdirectory(pth-mag)

PYTHONPATH=$1/build/src/ python $1/tests/unit/atomistic_anisotropy.py
PYTHONPATH=$1/build/src/ python $1/tests/unit/atomistic_dipole_dipole.py
PYTHONPATH=$1/build/src/ python $1/tests/unit/atomistic_exchange.py
PYTHONPATH=$1/build/src/ python $1/tests/unit/atomistic_dmi.py

#Planned:
#tests for all interactions
#tests for integration methods during refactoring
#"integration tests" for sp4, maybe sp5
