#!/bin/bash

# building tests
../../../scripts/build.sh .

# running all test binaries
for file in bin/*; do $file; done
