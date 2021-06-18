#!/bin/bash
# call into git toplevel, where compile_comands.json resides
cd $(git rev-parse --show-toplevel)

run-clang-tidy -checks='-*,clang-analyzer-*,performance-*' # -fix
# run-clang-tidy -checks='-*,modernize*' # -fix
# run-clang-tidy -checks='-*,modernize-use-override' -fix
