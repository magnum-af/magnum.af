#!/bin/bash -e
# $1 is supposed to be /path/to/magnum.af

# running all test binaries
for file in $1/tests/unit/cpp/bin/*; do $file; done
