#!/bin/bash
# $1 is supposed to be /path/to/magnum.af

# building tests
$1/scripts/bash/build.sh $1/tests/integration/cpp/

# running all test binaries
for file in $1/tests/integration/cpp/bin/*; do $file; done
