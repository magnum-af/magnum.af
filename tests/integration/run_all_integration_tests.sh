#!/bin/bash -e
#usage: ./runall /path/to/gitdirectory(magnum.af)

# cpp
$1/tests/integration/cpp/runall.sh $1

command -v pip3 >/dev/null 2>&1 && pip3_output="$(pip3 show arrayfire)" # Note: output of 'pip3 show' is empty if package not found
if [ -n "$pip3_output" ]; then
    python="python3"
else
    python="python"
fi

for filename in $1/tests/integration/python/*.py; do
    PYTHONPATH=$1/build/src/ $python $filename
done
