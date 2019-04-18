#!/bin/bash -e
#usage: ./runall /path/to/gitdirectory(magnum.af)

# cpp
./unit/run_all_unit_tests.sh $1
./integration/run_all_integration_tests.sh $1
