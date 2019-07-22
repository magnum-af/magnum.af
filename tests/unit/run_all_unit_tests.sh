#!/bin/bash -e
# script running all unit tests: i.e. the previously compiled (!) binaries in ./cpp and the python scripts in ./python

cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # call this scripts directory

# switch between python2 and python3
command -v pip3 >/dev/null 2>&1 && pip3_output="$(pip3 show arrayfire)" # Note: output of 'pip3 show' is empty if package not found
if [ -n "$pip3_output" ]; then
    python="python3"
else
    python="python"
fi

# running all python tests
for filename in ./python/*.py; do
    PYTHONPATH=../../build/src/ $python $filename
done
