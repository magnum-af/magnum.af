#!/bin/bash -xe
# building with clang static analyzer

# call into directory ../ w.r.t. this script
cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" && cd ..

mkdir -p build_analyze && cd build_analyze

# Note: $CXX is ignored, scan-build does use clang-compiler anyway
scan-build cmake ..
scan-build make -j
