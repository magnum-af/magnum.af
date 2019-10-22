#!/bin/bash -e
# script running all unit and integration tests: i.e. the previously compiled (!) binaries in ./cpp and the python scripts in ./python
cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # call this scripts directory
../test/runall.sh
../python/test/run_all_unit_tests.sh
../python/test/run_all_integration_tests.sh
