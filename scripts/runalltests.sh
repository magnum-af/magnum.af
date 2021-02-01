#!/bin/bash -e
# script running all tests: i.e. the previously compiled (!) binaries in ../bin/ and the python scripts in ../python/test/
cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # call this scripts directory

# run compiled C++ tests:
for testbin in ../bin/test*; do
    echo "running test '$testbin':"
    ./$testbin
    echo ""
done

# run python tests
for filename in ../python/test/*.py; do
    echo "running test '$filename'"
    PYTHONPATH=../build/python/ python3 $filename
done
