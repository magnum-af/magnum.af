#!/bin/bash
# $1 is supposed to be /path/to/magnum.af
set -e

# building tests
$1/scripts/bash/clean_build.sh $1/tests/integration/cpp/
