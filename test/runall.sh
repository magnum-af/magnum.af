#!/bin/bash -e
# script running all cpp tests: i.e. the previously compiled (!) binaries in ../bin/test* 
cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # call this scripts directory

for testbin in ../bin/test*; do
    echo "running test '$testbin':"
    ./$testbin
    echo ""
done
