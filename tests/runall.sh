#!/bin/bash
#usage: ./runall /path/to/gitdirectory(pth-mag)

#unit
for filename in $1/tests/unit/*.py; do
    PYTHONPATH=$1/build/src/ python $filename
done

#integration
for filename in $1/tests/integration/*.py; do
    PYTHONPATH=$1/build/src/ python $filename
done
