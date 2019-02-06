#!/bin/bash
#usage: ./runall /path/to/gitdirectory(magnum.af)

# cpp
$1/tests/unit/cpp/runall.sh $1
$1/tests/integration/cpp/runall.sh $1

# py unit
for filename in $1/tests/unit/*.py; do
    PYTHONPATH=$1/build/src/ python $filename
done

# py integration
for filename in $1/tests/integration/*.py; do
    PYTHONPATH=$1/build/src/ python $filename
done
