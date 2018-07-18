#!/bin/bash
# $1 is supposed to be /path/to/magnum.af

# building tests
$1/tests/integration/cpp/maketests.sh $1

# running all test binaries
$1/tests/integration/cpp/runall.sh $1
