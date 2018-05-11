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

#Planned:
#tests for all interactions
#tests for integration methods during refactoring
#"integration tests" for sp4, maybe sp5