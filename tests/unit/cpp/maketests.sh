#!/bin/bash -e
# $1 is supposed to be /path/to/magnum.af

# building tests
$1/scripts/bash/clean_build.sh $1/tests/unit/cpp/
