#!/bin/bash
#usage: ./runall /path/to/gitdirectory(magnum.af)

#unit
for filename in $1/tests/unit/*.py; do
    PYTHONPATH=$1/build/src/ python $filename
done

#integration
for filename in $1/tests/integration/*.py; do
    PYTHONPATH=$1/build/src/ python $filename
done
