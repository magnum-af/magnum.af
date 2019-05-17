#!/bin/bash -e
# builds and runns all cpp unit tests

cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # call this scripts directory

# building tests
./maketests.sh

# running all test binaries
./runall.sh
