#!/bin/bash -xe
# building with clang static analyzer

mkdir -p build_analyze && cd build_analyze
scan-build cmake ..
scan-build make -j
