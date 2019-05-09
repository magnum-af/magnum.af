#!/bin/bash -e
#usage: ./runall /path/to/gitdirectory(magnum.af)

# cpp
$1/tests/unit/run_all_unit_tests.sh $1
$1/tests/integration/run_all_integration_tests.sh $1
