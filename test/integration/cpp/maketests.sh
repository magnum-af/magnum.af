#!/bin/bash -e
cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # call this scripts directory

# building tests
../../../scripts/bash/clean_build.sh ../../../test/integration/cpp/
